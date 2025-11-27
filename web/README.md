# 8-band k·p Quantum Well Simulator - Web Interface

Modern web interface for configuring quantum well structures and visualizing bandstructure and g-factor calculations.

## Features

- **Material Selection**: Browse 40+ III-V semiconductor materials with searchable interface
- **Quantum Well Builder**: Visual layer-by-layer construction with real-time diagram
- **Calculation Settings**: Configure wave vectors, k-points, bands, and g-factor methods
- **Interactive Visualization**: Plotly.js-powered bandstructure and g-factor plots
- **Input File Generation**: One-click download of properly formatted `.example` files

## Quick Start

### 1. Start Local Server

```bash
cd /data/8bandkp-fdm-ai/web
python3 -m http.server 8000
```

### 2. Open in Browser

Navigate to: `http://localhost:8000`

### 3. Configure Your System

1. **Select Materials**: Choose from the material library (e.g., GaAs, AlGaAs, InAs)
2. **Build Structure**: Add/modify layers with start/end positions and offsets
3. **Set Parameters**: Configure wave vector, k-points, bands, etc.
4. **Generate Input**: Click "Generate Input" to download `.example` file

### 4. Run Simulation

```bash
cd /data/8bandkp-fdm-ai
./bandStructure < quantum_well.example
# or
./gfactorCalculation < quantum_well.example
```

### 5. Visualize Results

Upload the generated files in the visualization panel:
- `eigenvalues.dat` → Bandstructure plot
- `gfactor.dat` → G-factor bar chart

## Interface Layout

### Left Sidebar: Material Selection
- Searchable material library
- Material parameter preview (Eg, m*, Δ, γ₁, γ₂, γ₃)
- System type selector (Bulk/Quantum Well)

### Center Panel: Structure & Settings
- **Layer Builder**: Visual quantum well construction
- **Diagram**: Real-time structure visualization
- **Settings**: Wave vector, k-points, bands, g-factor method

### Right Panel: Visualization
- **Bandstructure**: Energy dispersion curves vs k-vector
- **G-Factors**: Component bar chart (gx, gy, gz)

## Material Database

Includes 40+ materials from `parameters.f90`:
- **III-V Binaries**: GaAs, InAs, AlAs, GaSb, AlSb, InSb, InP
- **AlGaAs Alloys**: Al15Ga85As, Al20Ga80As, Al30Ga70As
- **InGaAs Alloys**: In20Ga80As
- **InAsSb Alloys**: InAs10Sb90 through InAs90Sb10
- **Others**: Al63In37Sb

Each material displays:
- Band gap (Eg) in eV
- Effective mass (m*) in m₀
- Spin-orbit splitting (Δ) in eV
- Luttinger parameters (γ₁, γ₂, γ₃)

## Example Workflows

### Creating a GaAs/AlGaAs Quantum Well

1. Click "Reset" to start fresh
2. Default structure loads with:
   - Barrier: Al30Ga70As (-200Å to -50Å)
   - Well: GaAs (-50Å to 50Å)
   - Barrier: Al30Ga70As (50Å to 200Å)
3. Adjust parameters as needed
4. Generate and run simulation

### Calculating G-Factors

1. Build your structure
2. Set k-vector direction to "k0 (Γ-point)"
3. Set k-points to 0 or 1
4. Choose g-factor method:
   - **Analytical**: Works for bulk and QW
   - **Numerical**: Bulk systems only
5. Select band type (CB/VB) and index
6. Generate input and run `gfactorCalculation`

## Design Features

- **Dark Theme**: Professional navy background (#0f172a)
- **Glassmorphism**: Frosted glass effect on cards
- **Vibrant Accents**: Purple (#8b5cf6) and blue (#3b82f6) gradients
- **Smooth Animations**: Hover effects and transitions
- **Responsive Layout**: Adapts to different screen sizes
- **Material Colors**: Each semiconductor has a unique color

## Browser Compatibility

- Chrome/Edge (recommended)
- Firefox
- Safari
- Requires modern JavaScript (ES6+)

## Technical Stack

- **HTML5**: Semantic structure
- **CSS3**: Custom properties, Grid, Flexbox, animations
- **JavaScript**: ES6+ modules, classes
- **Plotly.js**: Interactive 2D graphing library
- **Google Fonts**: Inter font family

## Troubleshooting

**Blank plots after loading data:**
- Ensure file format matches expected structure
- Check browser console for errors

**Material list not showing:**
- Check browser JavaScript console
- Ensure `app.js` loaded correctly

**Download not working:**
- Check browser download settings
- Allow downloads from localhost

## Future Enhancements

- [ ] Backend integration for live simulation execution
- [ ] Wavefunction envelope plotting
- [ ] Batch simulation queue
- [ ] Export plots as PNG/SVG
- [ ] Material parameter editor
- [ ] Strain calculator
- [ ] Database of pre-configured examples

## License

Same as parent project (GPLv3)
