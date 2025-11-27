#!/usr/bin/env python3
"""
Flask backend server for 8-band k·p Quantum Well Simulator
Handles simulation execution and result retrieval
"""

from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
import subprocess
import os
import time
import uuid
import tempfile
from pathlib import Path

app = Flask(__name__, static_folder='.')
CORS(app)

# Base directory for the simulator
BASE_DIR = Path('/data/8bandkp-fdm-ai')
WEB_DIR = BASE_DIR / 'web'
RUN_DIR = BASE_DIR / 'run'

# Store running simulations
simulations = {}

@app.route('/')
def index():
    """Serve the main interface"""
    return send_from_directory('.', 'index.html')

@app.route('/<path:path>')
def serve_static(path):
    """Serve static files"""
    return send_from_directory('.', path)

@app.route('/api/run-simulation', methods=['POST'])
def run_simulation():
    """
    Run a simulation with the provided configuration
    
    Expected JSON:
    {
        "layers": [...],
        "settings": {...},
        "simulationType": "bandstructure" or "gfactor"
    }
    """
    try:
        data = request.json
        simulation_id = str(uuid.uuid4())[:8]
        
        # Create temporary directory for this simulation
        sim_dir = RUN_DIR / f'web_sim_{simulation_id}'
        sim_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate input file
        input_content = generate_input_file(data)
        input_file = sim_dir / 'input.cfg'
        input_file.write_text(input_content)
        
        # Determine which executable to run
        sim_type = data.get('simulationType', 'bandstructure')
        executable = BASE_DIR / ('gfactorCalculation' if sim_type == 'gfactor' else 'bandStructure')
        
        if not executable.exists():
            return jsonify({'error': f'Executable not found: {executable}'}), 500
        
        # Run simulation in background
        process = subprocess.Popen(
            [str(executable)],
            stdin=open(input_file),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=sim_dir,
            text=True
        )
        
        # Store simulation info
        simulations[simulation_id] = {
            'process': process,
            'dir': sim_dir,
            'type': sim_type,
            'status': 'running',
            'start_time': time.time()
        }
        
        return jsonify({
            'simulationId': simulation_id,
            'status': 'running',
            'message': 'Simulation started successfully'
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/simulation-status/<simulation_id>', methods=['GET'])
def get_simulation_status(simulation_id):
    """Check simulation status and retrieve results if complete"""
    
    if simulation_id not in simulations:
        return jsonify({'error': 'Simulation not found'}), 404
    
    sim_info = simulations[simulation_id]
    process = sim_info['process']
    sim_dir = sim_info['dir']
    
    # Check if process is still running
    poll_result = process.poll()
    
    if poll_result is None:
        # Still running
        elapsed = time.time() - sim_info['start_time']
        return jsonify({
            'status': 'running',
            'elapsed': round(elapsed, 1),
            'message': f'Simulation running for {elapsed:.1f}s...'
        })
    
    # Process completed
    if poll_result != 0:
        # Error occurred
        stderr = process.stderr.read()
        return jsonify({
            'status': 'error',
            'exitCode': poll_result,
            'error': stderr
        }), 500
    
    # Success - read results
    try:
        results = {}
        
        # Read eigenvalues
        eigenvalues_file = sim_dir / 'eigenvalues.dat'
        
        # Check if results are in outputs directory
        if not eigenvalues_file.exists():
            outputs_dir = sim_dir / 'outputs'
            if outputs_dir.exists():
                # Get latest directory
                subdirs = [d for d in outputs_dir.iterdir() if d.is_dir()]
                if subdirs:
                    latest_dir = max(subdirs, key=lambda p: p.stat().st_mtime)
                    eigenvalues_file = latest_dir / 'eigenvalues.dat'
                    
        if eigenvalues_file.exists():
            results['eigenvalues'] = eigenvalues_file.read_text()
        
        # Read g-factors if available
        gfactor_file = sim_dir / 'gfactor.dat'
        
        # Check if results are in outputs directory (if we found a latest_dir before, use it)
        if not gfactor_file.exists() and 'latest_dir' in locals():
            gfactor_file = latest_dir / 'gfactor.dat'
            
        if gfactor_file.exists():
            results['gfactor'] = gfactor_file.read_text()
        
        # Read stdout
        stdout = process.stdout.read()
        
        sim_info['status'] = 'complete'
        
        return jsonify({
            'status': 'complete',
            'results': results,
            'output': stdout,
            'message': 'Simulation completed successfully'
        })
        
    except Exception as e:
        return jsonify({'error': f'Error reading results: {str(e)}'}), 500

def generate_input_file(data):
    """Generate input file content from configuration"""
    
    layers = data.get('layers', [])
    settings = data.get('settings', {})
    
    content = '# Generated by 8-band k·p Simulator Web Interface\n'
    content += f'# Generated: {time.strftime("%Y-%m-%d %H:%M:%S")}\n\n'
    
    # Wave vector settings
    content += f"waveVector {settings.get('waveVector', 'kx')}\n"
    content += f"waveVectorMax {settings.get('waveVectorMax', 0.1)}\n"
    content += f"waveVectorStep {settings.get('waveVectorStep', 11)}\n"
    
    # System configuration
    confinement = settings.get('confinement', 1)
    content += f"confinement {confinement}\n"
    content += f"FDstep {settings.get('fdStep', 101)}\n"
    content += f"numLayers {len(layers)}\n"
    
    # Materials
    for i, layer in enumerate(layers):
        if confinement == 0:
            content += f"material{i+1} {layer['material']}\n"
        else:
            content += f"material{i+1} {layer['material']} {layer['start']} {layer['end']} {layer['offset']}\n"
    
    # Bands
    content += f"numcb {settings.get('numCB', 32)}\n"
    content += f"numvb {settings.get('numVB', 32)}\n"
    
    # External field
    content += f"ExternalField 0 EF\n"
    content += f"EFParams {settings.get('efParams', 0.0005)}\n"
    
    # G-factor settings
    gfactor_method = settings.get('gfactorMethod', 'none')
    if gfactor_method != 'none':
        content += f"\n# G-Factor Settings\n"
        content += f"gfactorMethod {gfactor_method}\n"
        content += f"gfactorBand {settings.get('gfactorBandType', 'cb')} {settings.get('gfactorBandIndex', 1)}\n"
    
    return content

if __name__ == '__main__':
    print(f"Starting server on port 8080...")
    print(f"Open http://localhost:8080 in your browser")
    print(f"BASE_DIR: {BASE_DIR}")
    print(f"WEB_DIR: {WEB_DIR}")
    print(f"RUN_DIR: {RUN_DIR}")
    
    # Ensure run directory exists
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    
    # Change to web directory
    os.chdir(WEB_DIR)
    
    # Run server
    app.run(host='0.0.0.0', port=8080, debug=True, threaded=True)
