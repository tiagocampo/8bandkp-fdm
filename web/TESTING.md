# How to Test the Web Interface

## Quick Test Instructions

### Step 1: Ensure the Web Server is Running
The web server should already be running at `http://localhost:8000`

### Step 2: Test Files Location
Test data files have been copied to the `web/` directory:
- `web/test_eigenvalues.dat` - Bandstructure data
- `web/test_gfactor.dat` - G-factor data

### Step 3: Upload and Visualize

1. **Open the interface**: `http://localhost:8000`

2. **Test Bandstructure Plot**:
   - Scroll to the right panel "Bandstructure" section
   - Click the "Load Data" button
   - Select `test_eigenvalues.dat` from the file picker
   - The plot should appear showing energy bands vs k-vector

3. **Test G-Factor Plot**:
   - Scroll to "G-Factors" section below
   - Click the "Load Data" button  
   - Select `test_gfactor.dat` from the file picker
   - Bar chart should appear with gx, gy, gz values

### Step 4: Generate Your Own Input File

1. Modify the quantum well structure in the center panel
2. Adjust calculation settings
3. Click "Generate Input" at the top
4. Save the `quantum_well.example` file

### Step 5: Run a Simulation

To get fresh data:
```bash
cd /data/8bandkp-fdm-ai

# Copy one of the existing examples
cp examples/qw_gaas_algaas.example input.cfg

# Run the simulation
./bandStructure < input.cfg

# Results will be in eigenvalues.dat
# Upload this file to the web interface
```

## Troubleshooting

**Plot doesn't appear after uploading:**
- Open browser console (F12) and check for errors
- Verify the data file format is correct
- Try refreshing the page and uploading again

**File picker doesn't open:**
- The file input might be hidden - look for the "Load Data" button
- Try clicking directly on the text "Load Data"

**Simulation fails to run:**
- Make sure you're piping input correctly: `./bandStructure < input_file.example`
- Check that the executables are built: `make all`
- Verify input file format matches the examples

## Current Known Issues

The g-factor data (test_gfactor.dat) shows all values as 4.0, which suggests it's from an incomplete run. For better testing:
```bash
# Run a proper g-factor calculation
cp examples/gfactor_qw_gaas_algaas_cb_analytical.example input.cfg
./gfactorCalculation < input.cfg
# Upload the resulting gfactor.dat
```
