# Troubleshooting Guide: 8bandkp-fdm-ai

**Comprehensive guide for resolving common issues and errors**

## Build Issues

### Issue 1: Fortran Compiler Not Found

**Error Message**:
```
gfortran: command not found
make: gfortran: No such file or directory
```

**Solution**:
```bash
# Ubuntu/Debian
sudo apt install gfortran

# CentOS/RHEL
sudo yum install gcc-gfortran

# Arch/Manjaro
sudo pacman -S gcc-fortran

# macOS
brew install gfortran
```

**Verification**:
```bash
gfortran --version
# Expected: GNU Fortran (GCC) 15.2.1+
```

### Issue 2: MKL Library Not Found

**Error Message**:
```
/usr/bin/ld: cannot find -lmkl_gf_lp64
/usr/bin/ld: cannot find -lmkl_intel_thread
```

**Solution**:
1. **Install Intel OneAPI Base Toolkit**:
   ```bash
   # Download from Intel website
   wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/...
   sudo sh ./l_BaseKit_p_2025.0.1_offline.sh
   ```

2. **Source MKL environment**:
   ```bash
   source /opt/intel/oneapi/setvars.sh
   ```

3. **Update Makefile paths** (if needed):
   ```makefile
   LDFLAGS=-I/usr/include -L/usr/lib64 -L/opt/intel/oneapi/mkl/2025.0/lib/intel64 \
           -I/opt/intel/oneapi/mkl/2025.0/include \
           -L/opt/intel/oneapi/compiler/2025.0/lib \
           -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_intel_thread \
           -lmkl_core -liomp5 -lpthread -lm -ldl -lfftw3
   ```

### Issue 3: FFTW3 Not Found

**Error Message**:
```
/usr/bin/ld: cannot find -lfftw3
```

**Solution**:
```bash
# Ubuntu/Debian
sudo apt install libfftw3-dev

# CentOS/RHEL
sudo yum install fftw3-devel

# Arch/Manjaro
sudo pacman -S fftw

# macOS
brew install fftw
```

**Verification**:
```bash
pkg-config --exists fftw3 && echo "FFTW3 found"
```

### Issue 4: OpenMP Library Not Found

**Error Message**:
```
/usr/bin/ld: cannot find -liomp5
```

**Solution**:
1. **Check OpenMP library location**:
   ```bash
   find /opt/intel -name "libiomp5.so" 2>/dev/null
   ```

2. **Update Makefile**:
   ```makefile
   LDFLAGS=... -L/opt/intel/oneapi/compiler/2025.0/lib -liomp5 ...
   ```

3. **Alternative: Use GNU OpenMP**:
   ```makefile
   LDFLAGS=... -lgomp ...
   ```

### Issue 5: Undefined MKL Sparse Functions

**Error Message**:
```
undefined reference to `MKL_SPARSE_Z_MV'
undefined reference to `MKL_SPARSE_Z_CREATE_COO'
```

**Solution**:
1. **Ensure MKL sparse libraries are linked**:
   ```makefile
   LDFLAGS=... -lmkl_gf_lp64 -lmkl_intel_thread -lmkl_core ...
   ```

2. **Check MKL installation**:
   ```bash
   ls /opt/intel/oneapi/mkl/2025.0/lib/intel64/libmkl_*.so
   ```

## Runtime Issues

### Issue 6: Input File Not Found

**Error Message**:
```
STOP: Cannot open input file
```

**Solution**:
1. **Check file exists**:
   ```bash
   ls -la examples/bulk_InAs60Sb40.example
   ```

2. **Check file permissions**:
   ```bash
   chmod 644 examples/bulk_InAs60Sb40.example
   ```

3. **Verify file format**:
   ```bash
   head -5 examples/bulk_InAs60Sb40.example
   ```

### Issue 7: G-Factor Calculation Requires k=0

**Error Message**:
```
STOP g-factor calculation requires only k=0
```

**Solution**:
1. **Use correct input format**:
   ```
   waveVector: k0
   waveVectorStep: 0
   ```

2. **Check example file**:
   ```bash
   cat examples/gfactor_quantum_well.example
   ```

### Issue 8: Memory Issues

**Error Message**:
```
Segmentation fault
Out of memory
```

**Solution**:
1. **Reduce problem size**:
   ```
   FDstep: 51        # Instead of 101
   numcb: 16         # Instead of 32
   numvb: 16         # Instead of 32
   ```

2. **Check system memory**:
   ```bash
   free -h
   ```

3. **Use swap if available**:
   ```bash
   sudo swapon -a
   ```

### Issue 9: Numerical Convergence Issues

**Error Message**:
```
STOP: Matrix diagonalization failed
STOP: Eigenvalue calculation did not converge
```

**Solution**:
1. **Check input parameters**:
   - Ensure material parameters are reasonable
   - Verify band offsets are correct
   - Check discretization is sufficient

2. **Adjust numerical parameters**:
   ```
   FDstep: 201       # Increase discretization
   ```

3. **Check for unphysical values**:
   ```bash
   grep -i "nan\|inf" eigenvalues.dat
   ```

## Output Issues

### Issue 10: No Output Files Generated

**Symptoms**:
- Calculation runs but no files created
- Missing eigenvalues.dat

**Solution**:
1. **Check outputs directory**:
   ```bash
   ls -la outputs/*/*
   ```

2. **Check file permissions**:
   ```bash
   ls -la
   chmod 755 .
   ```

3. **Run with verbose output**:
   ```bash
   ./bandStructure examples/bulk_InAs60Sb40.example --out outputs/$(date +%Y%m%d-%H%M%S)/bulk 2>&1 | tee output.log
   ```

### Issue 11: Unphysical Results

**Symptoms**:
- Energy values outside reasonable range
- Band gap too large or too small
- No quantum well effects

**Solution**:
1. **Check material parameters**:
   ```bash
   grep -A 10 "Material:" output.log
   ```

2. **Verify input file format**:
   ```bash
   cat examples/bulk_InAs60Sb40.example
   ```

3. **Compare with validation results**:
   ```bash
   ./scripts/validate_results.sh
   ```

### Issue 12: Plotting Issues

**Error Message**:
```
gnuplot: command not found
Cannot open data file
```

**Solution**:
1. **Install gnuplot**:
   ```bash
   sudo apt install gnuplot  # Ubuntu/Debian
   sudo pacman -S gnuplot    # Arch/Manjaro
   brew install gnuplot      # macOS
   ```

2. **Check data file exists**:
   ```bash
   ls -la eigenvalues.dat
   ```

3. **Test gnuplot**:
   ```bash
   echo "plot sin(x)" | gnuplot
   ```

## Performance Issues

### Issue 13: Slow Build

**Symptoms**:
- Build takes > 10 minutes
- Compilation hangs

**Solution**:
1. **Use parallel compilation**:
   ```bash
   make -j4 all
   ```

2. **Check system resources**:
   ```bash
   top
   htop
   ```

3. **Disable optimizations temporarily**:
   ```makefile
   FFLAGS=-O0 -g
   ```

### Issue 14: Slow Calculations

**Symptoms**:
- Calculation takes > 5 minutes
- High CPU usage

**Solution**:
1. **Reduce problem size**:
   ```
   FDstep: 51
   waveVectorStep: 5
   ```

2. **Check system performance**:
   ```bash
   iostat -x 1
   ```

3. **Use optimized libraries**:
   - Ensure MKL is properly linked
   - Check OpenMP threading

## System-Specific Issues

### Issue 15: Ubuntu/Debian Specific

**Error**: Package conflicts
**Solution**:
```bash
sudo apt update
sudo apt upgrade
sudo apt install gfortran gcc g++ gnuplot make
```

### Issue 16: CentOS/RHEL Specific

**Error**: Missing development tools
**Solution**:
```bash
sudo yum groupinstall "Development Tools"
sudo yum install gcc-gfortran gnuplot make
```

### Issue 17: Arch/Manjaro Specific

**Error**: AUR package issues
**Solution**:
```bash
sudo pacman -S base-devel
yay -S intel-oneapi-basekit
```

### Issue 18: macOS Specific

**Error**: Xcode command line tools
**Solution**:
```bash
xcode-select --install
brew install gfortran gnuplot
```

## Debugging Tools

### Verbose Output
```bash
# Build with verbose output
make all V=1

# Run with debug output
./bandStructure input.example --out outputs/$(date +%Y%m%d-%H%M%S)/debug 2>&1 | tee debug.log
```

### Memory Debugging
```bash
# Check memory usage
valgrind --tool=memcheck ./bandStructure input.example --out outputs/valgrind-run

# Monitor memory
watch -n 1 'ps aux | grep bandStructure'
```

### Performance Profiling
```bash
# Profile with gprof
gfortran -pg -O2 -c src/*.f90
make all
./bandStructure input.example --out outputs/profile-run
gprof bandStructure gmon.out > profile.txt
```

## Getting Help

### 1. Check Documentation
- [Build Guide](BUILD.md)
- [Dependencies](DEPENDENCIES.md)
- [Validation](VALIDATION.md)
- [Plotting](PLOTTING.md)

### 2. Run Diagnostic Scripts
```bash
# Check system requirements
./scripts/verify_build.sh

# Validate results
./scripts/validate_results.sh
```

### 3. Create Minimal Test Case
```bash
# Create simple test
cat > test_simple.example << EOF
waveVector: k0
waveVectorMax: 0.01
waveVectorStep: 1
confinement: 0
FDstep: 1
numLayers: 1
material3: InAs60Sb40
numcb: 2
numvb: 6
ExternalField: 0 EF
EFParams: 0.0005
EOF

# Test with minimal case
./bandStructure < test_simple.example
```

### 4. Report Issues
When reporting issues, include:
- Operating system and version
- Compiler version (`gfortran --version`)
- Complete error message
- Input file used
- System specifications
- Steps to reproduce

## Prevention

### Best Practices
1. **Always use example files** as starting points
2. **Validate results** after each calculation
3. **Keep system updated** with latest packages
4. **Use version control** for input files
5. **Document changes** to input parameters

### Regular Maintenance
1. **Clean build artifacts** regularly:
   ```bash
   make clean_all
   ```

2. **Update dependencies** periodically
3. **Run validation tests** before important calculations
4. **Backup important results**

## Success Indicators

You know everything is working correctly when:
- ✅ Build completes without errors
- ✅ Example calculations run successfully
- ✅ Results match validation criteria
- ✅ Plots generate correctly
- ✅ No error messages in output
- ✅ Physical results are reasonable
