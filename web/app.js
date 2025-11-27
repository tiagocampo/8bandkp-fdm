// Material Database (extracted from parameters.f90)
const MATERIALS = {
    'GaAs': { Eg: 1.519, meff: 0.067, deltaSO: 0.341, gamma1: 6.98, gamma2: 2.06, gamma3: 2.93, color: '#3b82f6' },
    'Al20Ga80As': { Eg: 1.835, meff: 0.0836, deltaSO: 0.3288, gamma1: 6.336, gamma2: 1.812, gamma3: 2.628, color: '#8b5cf6' },
    'Al15Ga85As': { Eg: 1.756, meff: 0.0795, deltaSO: 0.3319, gamma1: 6.497, gamma2: 1.874, gamma3: 2.7035, color: '#7c3aed' },
    'Al30Ga70As': { Eg: 1.977, meff: 0.093, deltaSO: 0.353, gamma1: 6.107, gamma2: 1.773, gamma3: 2.543, color: '#9333ea' },
    'In20Ga80As': { Eg: 1.299, meff: 0.0607, deltaSO: 0.351, gamma1: 9.604, gamma2: 3.48, gamma3: 4.174, color: '#06b6d4' },
    'InAs': { Eg: 0.417, meff: 0.026, deltaSO: 0.39, gamma1: 20.0, gamma2: 8.5, gamma3: 9.2, color: '#10b981' },
    'AlAs': { Eg: 3.099, meff: 0.15, deltaSO: 0.28, gamma1: 3.76, gamma2: 0.82, gamma3: 1.42, color: '#a855f7' },
    'InP': { Eg: 1.4236, meff: 0.0795, deltaSO: 0.108, gamma1: 5.08, gamma2: 1.60, gamma3: 2.10, color: '#14b8a6' },
    'GaSb': { Eg: 0.812, meff: 0.039, deltaSO: 0.76, gamma1: 13.4, gamma2: 4.7, gamma3: 6.0, color: '#f59e0b' },
    'AlSb': { Eg: 2.386, meff: 0.14, deltaSO: 0.676, gamma1: 5.18, gamma2: 1.19, gamma3: 1.97, color: '#d946ef' },
    'InSb': { Eg: 0.235, meff: 0.0135, deltaSO: 0.81, gamma1: 34.8, gamma2: 15.5, gamma3: 16.5, color: '#ef4444' },
    'InAs10Sb90': { Eg: 0.1749, meff: 0.0116, deltaSO: 0.66, gamma1: 33.32, gamma2: 14.8, gamma3: 15.77, color: '#ec4899' },
    'InAs20Sb80': { Eg: 0.1322, meff: 0.0104, deltaSO: 0.534, gamma1: 31.84, gamma2: 14.1, gamma3: 15.04, color: '#f43f5e' },
    'InAs30Sb70': { Eg: 0.1069, meff: 0.0099, deltaSO: 0.432, gamma1: 30.36, gamma2: 13.4, gamma3: 14.31, color: '#fb7185' },
    'InAs40Sb60': { Eg: 0.099, meff: 0.0101, deltaSO: 0.354, gamma1: 28.88, gamma2: 12.7, gamma3: 13.58, color: '#fda4af' },
    'InAs50Sb50': { Eg: 0.1085, meff: 0.011, deltaSO: 0.3, gamma1: 27.4, gamma2: 12.0, gamma3: 12.85, color: '#fb923c' },
    'InAs60Sb40': { Eg: 0.1354, meff: 0.0126, deltaSO: 0.27, gamma1: 25.92, gamma2: 11.30, gamma3: 12.12, color: '#fdba74' },
    'InAs70Sb30': { Eg: 0.1797, meff: 0.0149, deltaSO: 0.264, gamma1: 24.44, gamma2: 10.60, gamma3: 11.39, color: '#fbbf24' },
    'InAs80Sb20': { Eg: 0.2414, meff: 0.0179, deltaSO: 0.282, gamma1: 22.96, gamma2: 9.90, gamma3: 10.66, color: '#fcd34d' },
    'InAs90Sb10': { Eg: 0.3205, meff: 0.0216, deltaSO: 0.324, gamma1: 21.48, gamma2: 9.20, gamma3: 9.93, color: '#fde047' },
    'Al63In37Sb': { Eg: 0.9306, meff: 0.0137, deltaSO: 0.7021, gamma1: 21.7616, gamma2: 9.0398, gamma3: 10.0529, color: '#c084fc' },
};

// Backend API base URL
const API_BASE = 'http://localhost:8080/api';

// Application State
const app = {
    selectedMaterial: null,
    layers: [],
    confinement: 1,
    currentSimulation: null,
    pollInterval: null,
    settings: {
        waveVector: 'kx',
        waveVectorMax: 0.1,
        waveVectorStep: 11,
        fdStep: 101,
        numCB: 32,
        numVB: 32,
        efParams: 0.0005,
        gfactorMethod: 'none',
        gfactorBandType: 'cb',
        gfactorBandIndex: 1
    },

    init() {
        this.renderMaterialList();
        this.addDefaultLayers();
        this.updateDiagram();
        this.initPlots();
    },

    renderMaterialList() {
        const listEl = document.getElementById('materialList');
        const materials = Object.keys(MATERIALS).sort();
        
        listEl.innerHTML = materials.map(name => `
            <div class="material-item" onclick="app.selectMaterial('${name}')">
                <div class="material-dot" style="background: ${MATERIALS[name].color}"></div>
                <div class="material-name">${name}</div>
            </div>
        `).join('');
    },

    filterMaterials(query) {
        const listEl = document.getElementById('materialList');
        const materials = Object.keys(MATERIALS).sort();
        const filtered = materials.filter(name => 
            name.toLowerCase().includes(query.toLowerCase())
        );
        
        listEl.innerHTML = filtered.map(name => `
            <div class="material-item" onclick="app.selectMaterial('${name}')">
                <div class="material-dot" style="background: ${MATERIALS[name].color}"></div>
                <div class="material-name">${name}</div>
            </div>
        `).join('');
    },

    selectMaterial(name) {
        this.selectedMaterial = name;
        
        // Update UI
        document.querySelectorAll('.material-item').forEach(el => {
            el.classList.remove('selected');
        });
        event.currentTarget.classList.add('selected');
        
        // Show material info
        const mat = MATERIALS[name];
        const infoEl = document.getElementById('materialInfo');
        infoEl.innerHTML = `
            <div class="info-row">
                <span class="info-label">Band Gap (Eg)</span>
                <span class="info-value">${mat.Eg.toFixed(3)} eV</span>
            </div>
            <div class="info-row">
                <span class="info-label">Eff. Mass (m*)</span>
                <span class="info-value">${mat.meff.toFixed(4)} m₀</span>
            </div>
            <div class="info-row">
                <span class="info-label">Spin-Orbit (Δ)</span>
                <span class="info-value">${mat.deltaSO.toFixed(3)} eV</span>
            </div>
            <div class="info-row">
                <span class="info-label">γ₁</span>
                <span class="info-value">${mat.gamma1.toFixed(2)}</span>
            </div>
            <div class="info-row">
                <span class="info-label">γ₂</span>
                <span class="info-value">${mat.gamma2.toFixed(2)}</span>
            </div>
            <div class="info-row">
                <span class="info-label">γ₃</span>
                <span class="info-value">${mat.gamma3.toFixed(2)}</span>
            </div>
        `;
    },

    addDefaultLayers() {
        this.layers = [
            { material: 'Al30Ga70As', start: -200, end: 200, offset: 0 },
            { material: 'GaAs', start: -50, end: 50, offset: -0.258 }
        ];
        this.renderLayers();
    },

    addLayer() {
        this.layers.push({
            material: 'GaAs',
            start: 0,
            end: 100,
            offset: 0
        });
        this.renderLayers();
        this.updateDiagram();
    },

    removeLayer(index) {
        this.layers.splice(index, 1);
        this.renderLayers();
        this.updateDiagram();
    },

    updateLayer(index, field, value) {
        this.layers[index][field] = field === 'material' ? value : parseFloat(value);
        this.updateDiagram();
    },

    renderLayers() {
        const builderEl = document.getElementById('layerBuilder');
        const materials = Object.keys(MATERIALS).sort();
        
        builderEl.innerHTML = this.layers.map((layer, i) => `
            <div class="layer-row">
                <div class="layer-color" style="background: ${MATERIALS[layer.material].color}"></div>
                <select onchange="app.updateLayer(${i}, 'material', this.value)">
                    ${materials.map(m => `
                        <option value="${m}" ${m === layer.material ? 'selected' : ''}>${m}</option>
                    `).join('')}
                </select>
                <input type="number" value="${layer.start}" 
                       onchange="app.updateLayer(${i}, 'start', this.value)"
                       placeholder="Start (Å)">
                <input type="number" value="${layer.end}"
                       onchange="app.updateLayer(${i}, 'end', this.value)"
                       placeholder="End (Å)">
                <input type="number" value="${layer.offset}" step="0.001"
                       onchange="app.updateLayer(${i}, 'offset', this.value)"
                       placeholder="Offset">
                <button onclick="app.removeLayer(${i})">
                    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor">
                        <line x1="18" y1="6" x2="6" y2="18"></line>
                        <line x1="6" y1="6" x2="18" y2="18"></line>
                    </svg>
                </button>
            </div>
        `).join('');
    },

    updateDiagram() {
        const svg = document.getElementById('layerSvg');
        if (!this.layers.length) {
            svg.innerHTML = '<text x="400" y="150" text-anchor="middle" fill="#94a3b8">Add layers to visualize structure</text>';
            return;
        }

        // Calculate bounds
        const allPositions = this.layers.flatMap(l => [l.start, l.end]);
        const minPos = Math.min(...allPositions);
        const maxPos = Math.max(...allPositions);
        const range = maxPos - minPos;
        
        const width = 800;
        const height = 300;
        const padding = 50;
        const barHeight = 80;
        const barY = (height - barHeight) / 2;

        // Scale function
        const scale = (pos) => padding + ((pos - minPos) / range) * (width - 2 * padding);

        // Draw layers
        let layersHTML = '';
        this.layers.forEach((layer, i) => {
            const x = scale(layer.start);
            const w = scale(layer.end) - x;
            const color = MATERIALS[layer.material].color;
            
            layersHTML += `
                <rect x="${x}" y="${barY}" width="${w}" height="${barHeight}" 
                      fill="${color}" opacity="0.7" stroke="${color}" stroke-width="2" rx="4"/>
                <text x="${x + w/2}" y="${barY + barHeight/2}" 
                      text-anchor="middle" dominant-baseline="middle" 
                      fill="white" font-size="12" font-weight="600">
                    ${layer.material}
                </text>
                <text x="${x + w/2}" y="${barY + barHeight + 20}" 
                      text-anchor="middle" fill="#cbd5e1" font-size="11">
                    ${layer.start}Å to ${layer.end}Å
                </text>
            `;
        });

        // Draw axis
        const axisY = barY + barHeight + 50;
        layersHTML += `
            <line x1="${padding}" y1="${axisY}" x2="${width - padding}" y2="${axisY}" 
                  stroke="#475569" stroke-width="2"/>
            <text x="${width/2}" y="${axisY + 30}" text-anchor="middle" 
                  fill="#94a3b8" font-size="12">Position (Angstroms)</text>
        `;

        svg.innerHTML = layersHTML;
    },

    updateConfinement(value) {
        this.confinement = parseInt(value);
        if (this.confinement === 0) {
            // Bulk system
            this.layers = [{ material: 'GaAs', start: 0, end: 0, offset: 0 }];
            this.renderLayers();
            this.updateDiagram();
        }
    },

    updateSettings() {
        this.settings.waveVector = document.getElementById('waveVector').value;
        this.settings.waveVectorMax = parseFloat(document.getElementById('waveVectorMax').value);
        this.settings.waveVectorStep = parseInt(document.getElementById('waveVectorStep').value);
        this.settings.fdStep = parseInt(document.getElementById('fdStep').value);
        this.settings.numCB = parseInt(document.getElementById('numCB').value);
        this.settings.numVB = parseInt(document.getElementById('numVB').value);
        this.settings.efParams = parseFloat(document.getElementById('efParams').value);
        this.settings.gfactorMethod = document.getElementById('gfactorMethod').value;
        this.settings.gfactorBandType = document.getElementById('gfactorBandType')?.value || 'cb';
        this.settings.gfactorBandIndex = parseInt(document.getElementById('gfactorBandIndex')?.value || 1);

        // Show/hide g-factor settings
        const gfactorSettings = document.getElementById('gfactorSettings');
        if (this.settings.gfactorMethod !== 'none') {
            gfactorSettings.style.display = 'grid';
        } else {
            gfactorSettings.style.display = 'none';
        }
    },

    generateInputFile() {
        let content = '# Generated by 8-band k·p Simulator Web Interface\n';
        content += `# Generated: ${new Date().toISOString()}\n\n`;
        
        // Wave vector settings
        content += `waveVector ${this.settings.waveVector}\n`;
        content += `waveVectorMax ${this.settings.waveVectorMax}\n`;
        content += `waveVectorStep ${this.settings.waveVectorStep}\n`;
        
        // System configuration
        content += `confinement ${this.confinement}\n`;
        content += `FDstep ${this.settings.fdStep}\n`;
        content += `numLayers ${this.layers.length}\n`;
        
        // Materials
        this.layers.forEach((layer, i) => {
            if (this.confinement === 0) {
                content += `material${i+1} ${layer.material}\n`;
            } else {
                content += `material${i+1} ${layer.material} ${layer.start} ${layer.end} ${layer.offset}\n`;
            }
        });
        
        // Bands
        content += `numcb ${this.settings.numCB}\n`;
        content += `numvb ${this.settings.numVB}\n`;
        
        // External field
        content += `ExternalField 0 EF\n`;
        content += `EFParams ${this.settings.efParams}\n`;
        
        // G-factor settings
        if (this.settings.gfactorMethod !== 'none') {
            content += `\n# G-Factor Settings\n`;
            content += `gfactorMethod ${this.settings.gfactorMethod}\n`;
            content += `gfactorBand ${this.settings.gfactorBandType} ${this.settings.gfactorBandIndex}\n`;
        }
        
        // Download file
        const blob = new Blob([content], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'quantum_well.example';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
        
        // Show instructions
        alert(`Input file generated!\n\nTo run the simulation:\n1. Copy the file to your simulator directory\n2. Run: ./bandStructure < quantum_well.example\n   or: ./gfactorCalculation < quantum_well.example\n3. Upload the results using the "Load Data" buttons`);
    },

    initPlots() {
        // Initialize empty plots
        const layout = {
            paper_bgcolor: '#0f172a',
            plot_bgcolor: '#1e293b',
            font: { color: '#cbd5e1' },
            margin: { t: 40, r: 20, b: 40, l: 50 },
            xaxis: { gridcolor: '#334155' },
            yaxis: { gridcolor: '#334155' }
        };

        Plotly.newPlot('bandstructurePlot', [], {
            ...layout,
            title: 'Upload eigenvalues.dat to visualize',
            xaxis: { ...layout.xaxis, title: 'k-vector' },
            yaxis: { ...layout.yaxis, title: 'Energy (eV)' }
        }, { responsive: true });

        Plotly.newPlot('gfactorPlot', [], {
            ...layout,
            title: 'Upload gfactor.dat to visualize',
            xaxis: { ...layout.xaxis, title: 'Component' },
            yaxis: { ...layout.yaxis, title: 'g-factor' }
        }, { responsive: true });
    },

    loadEigenvalues(file) {
        const reader = new FileReader();
        reader.onload = (e) => {
            const text = e.target.result;
            const lines = text.trim().split('\n').filter(l => !l.startsWith('#'));
            
            const data = lines.map(line => {
                const values = line.trim().split(/\s+/).map(parseFloat);
                return values;
            });

            if (data.length === 0) return;

            const numBands = data[0].length - 1;
            const traces = [];
            const colors = ['#3b82f6', '#8b5cf6', '#10b981', '#f59e0b', '#ef4444', '#ec4899', '#06b6d4', '#a855f7'];

            for (let band = 0; band < numBands; band++) {
                traces.push({
                    x: data.map(row => row[0]),
                    y: data.map(row => row[band + 1]),
                    mode: 'lines',
                    name: `Band ${band + 1}`,
                    line: {
                        color: colors[band % colors.length],
                        width: 2
                    }
                });
            }

            const layout = {
                title: 'Band Structure',
                paper_bgcolor: '#0f172a',
                plot_bgcolor: '#1e293b',
                font: { color: '#cbd5e1', family: 'Inter' },
                margin: { t: 40, r: 20, b: 50, l: 60 },
                xaxis: {
                    title: 'k-vector',
                    gridcolor: '#334155',
                    color: '#cbd5e1'
                },
                yaxis: {
                    title: 'Energy (eV)',
                    gridcolor: '#334155',
                    color: '#cbd5e1'
                },
                showlegend: true,
                legend: {
                    x: 1.02,
                    y: 1,
                    bgcolor: 'rgba(30, 41, 59, 0.8)'
                }
            };

            Plotly.newPlot('bandstructurePlot', traces, layout, { responsive: true });
        };
        reader.readAsText(file);
    },

    loadGfactor(file) {
        const reader = new FileReader();
        reader.onload = (e) => {
            const text = e.target.result;
            const lines = text.trim().split('\n').filter(l => !l.startsWith('#'));
            
            if (lines.length === 0) return;
            
            const values = lines[0].trim().split(/\s+/).map(parseFloat);
            if (values.length < 3) return;

            const [gx, gy, gz] = values;

            const trace = {
                x: ['g_x', 'g_y', 'g_z'],
                y: [gx, gy, gz],
                type: 'bar',
                marker: {
                    color: ['#8b5cf6', '#3b82f6', '#10b981'],
                    line: {
                        color: '#1e293b',
                        width: 2
                    }
                },
                text: [gx.toFixed(3), gy.toFixed(3), gz.toFixed(3)],
                textposition: 'outside',
                textfont: { color: '#cbd5e1' }
            };

            const layout = {
                title: 'G-Factor Components',
                paper_bgcolor: '#0f172a',
                plot_bgcolor: '#1e293b',
                font: { color: '#cbd5e1', family: 'Inter' },
                margin: { t: 40, r: 20, b: 50, l: 60 },
                xaxis: {
                    title: 'Component',
                    gridcolor: '#334155',
                    color: '#cbd5e1'
                },
                yaxis: {
                    title: 'g-factor',
                    gridcolor: '#334155',
                    color: '#cbd5e1'
                },
                showlegend: false
            };

            Plotly.newPlot('gfactorPlot', [trace], layout, { responsive: true });
        };
        reader.readAsText(file);
    },

    async runSimulation() {
        try {            // Update settings first
            this.updateSettings();
            
            // Determine simulation type
            const simulationType = this.settings.gfactorMethod !== 'none' ? 'gfactor' : 'bandstructure';
            
            // Show loading state
            this.showSimulationStatus('Starting simulation...', 'running');
            
            // Send request to backend
            const response = await fetch(`${API_BASE}/run-simulation`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    layers: this.layers,
                    settings: { ...this.settings, confinement: this.confinement },
                    simulationType: simulationType
                })
            });
            
            if (!response.ok) {
                throw new Error(`Server error: ${response.statusText}`);
            }
            
            const data = await response.json();
            this.currentSimulation = data.simulationId;
            
            // Start polling for results
            this.pollSimulation();
            
        } catch (error) {
            console.error('Simulation error:', error);
            this.showSimulationStatus(`Error: ${error.message}`, 'error');
        }
    },

    pollSimulation() {
        // Clear any existing poll
        if (this.pollInterval) {
            clearInterval(this.pollInterval);
        }
        
        // Poll every 2 seconds
        this.pollInterval = setInterval(async () => {
            try {
                const response = await fetch(`${API_BASE}/simulation-status/${this.currentSimulation}`);
                const data = await response.json();
                
                if (data.status === 'running') {
                    this.showSimulationStatus(`Running... (${data.elapsed}s)`, 'running');
                } else if (data.status === 'complete') {
                    clearInterval(this.pollInterval);
                    this.handleSimulationComplete(data.results);
                } else if (data.status === 'error') {
                    clearInterval(this.pollInterval);
                    this.showSimulationStatus(`Error: ${data.error}`, 'error');
                }
            } catch (error) {
                clearInterval(this.pollInterval);
                this.showSimulationStatus(`Polling error: ${error.message}`, 'error');
            }
        }, 2000);
    },

    handleSimulationComplete(results) {
        this.showSimulationStatus('Simulation complete! Plotting results...', 'success');
        
        // Auto-plot eigenvalues if available
        if (results.eigenvalues) {
            this.plotEigenvaluesFromText(results.eigenvalues);
        }
        
        // Auto-plot g-factors if available
        if (results.gfactor) {
            this.plotGfactorFromText(results.gfactor);
        }
        
        setTimeout(() => {
            this.showSimulationStatus('', 'idle');
        }, 3000);
    },

    plotEigenvaluesFromText(text) {
        const lines = text.trim().split('\n').filter(l => !l.startsWith('#') && l.trim().length > 0);
        
        const data = lines.map(line => {
            const values = line.trim().split(/\s+/).map(parseFloat);
            return values;
        });

        if (data.length === 0) return;

        const numBands = data[0].length - 1;
        const traces = [];
        const colors = ['#3b82f6', '#8b5cf6', '#10b981', '#f59e0b', '#ef4444', '#ec4899', '#06b6d4', '#a855f7'];

        for (let band = 0; band < numBands; band++) {
            traces.push({
                x: data.map(row => row[0]),
                y: data.map(row => row[band + 1]),
                mode: 'lines',
                name: `Band ${band + 1}`,
                line: {
                    color: colors[band % colors.length],
                    width: 2
                }
            });
        }

        const layout = {
            title: 'Band Structure',
            paper_bgcolor: '#0f172a',
            plot_bgcolor: '#1e293b',
            font: { color: '#cbd5e1', family: 'Inter' },
            margin: { t: 40, r: 20, b: 50, l: 60 },
            xaxis: {
                title: 'k-vector',
                gridcolor: '#334155',
                color: '#cbd5e1'
            },
            yaxis: {
                title: 'Energy (eV)',
                gridcolor: '#334155',
                color: '#cbd5e1'
            },
            showlegend: true,
            legend: {
                x: 1.02,
                y: 1,
                bgcolor: 'rgba(30, 41, 59, 0.8)'
            }
        };

        Plotly.newPlot('bandstructurePlot', traces, layout, { responsive: true });
    },

    plotGfactorFromText(text) {
        const lines = text.trim().split('\n').filter(l => !l.startsWith('#') && l.trim().length > 0);
        
        if (lines.length === 0) return;
        
        const values = lines[0].trim().split(/\s+/).map(parseFloat);
        if (values.length < 3) return;

        const [gx, gy, gz] = values;

        const trace = {
            x: ['g_x', 'g_y', 'g_z'],
            y: [gx, gy, gz],
            type: 'bar',
            marker: {
                color: ['#8b5cf6', '#3b82f6', '#10b981'],
                line: {
                    color: '#1e293b',
                    width: 2
                }
            },
            text: [gx.toFixed(3), gy.toFixed(3), gz.toFixed(3)],
            textposition: 'outside',
            textfont: { color: '#cbd5e1' }
        };

        const layout = {
            title: 'G-Factor Components',
            paper_bgcolor: '#0f172a',
            plot_bgcolor: '#1e293b',
            font: { color: '#cbd5e1', family: 'Inter' },
            margin: { t: 40, r: 20, b: 50, l: 60 },
            xaxis: {
                title: 'Component',
                gridcolor: '#334155',
                color: '#cbd5e1'
            },
            yaxis: {
                title: 'g-factor',
                gridcolor: '#334155',
                color: '#cbd5e1'
            },
            showlegend: false
        };

        Plotly.newPlot('gfactorPlot', [trace], layout, { responsive: true });
    },

    showSimulationStatus(message, status) {
        // Update button text/state
        const runBtn = document.getElementById('runSimButton');
        if (runBtn) {
            if (status === 'running') {
                runBtn.disabled = true;
                runBtn.innerHTML = `
                    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" class="spinning">
                        <path d="M12 2v4M12 18v4M4.93 4.93l2.83 2.83M16.24 16.24l2.83 2.83M2 12h4M18 12h4M4.93 19.07l2.83-2.83M16.24 7.76l2.83-2.83"/>
                    </svg>
                    ${message}
                `;
            } else if (status === 'success') {
                runBtn.disabled = false;
                runBtn.innerHTML = `
                    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor">
                        <polyline points="20 6 9 17 4 12"/>
                    </svg>
                    ${message}
                `;
            } else if (status === 'error') {
                runBtn.disabled = false;
                runBtn.innerHTML = `
                    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor">
                        <circle cx="12" cy="12" r="10"/>
                        <line x1="15" y1="9" x2="9" y2="15"/>
                        <line x1="9" y1="9" x2="15" y2="15"/>
                    </svg>
                    ${message}
                `;
            } else {
                runBtn.disabled = false;
                runBtn.innerHTML = `
                    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor">
                        <polygon points="5 3 19 12 5 21 5 3"/>
                    </svg>
                    Run Simulation
                `;
            }
        }
    },

    resetAll() {
        if (confirm('Reset all settings to defaults?')) {
            document.getElementById('waveVector').value = 'kx';
            document.getElementById('waveVectorMax').value = 0.1;
            document.getElementById('waveVectorStep').value = 11;
            document.getElementById('fdStep').value = 101;
            document.getElementById('numCB').value = 32;
            document.getElementById('numVB').value = 32;
            document.getElementById('efParams').value = 0.0005;
            document.getElementById('gfactorMethod').value = 'none';
            
            this.addDefaultLayers();
            this.updateSettings();
            this.initPlots();
        }
    }
};

// Initialize on page load
document.addEventListener('DOMContentLoaded', () => {
    app.init();
});
