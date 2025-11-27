
text = """#k E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15 E16 
    0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
   0.100000E-01    0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
"""

lines = [l for l in text.strip().split('\n') if not l.startswith('#') and len(l.strip()) > 0]
print(f"Lines: {len(lines)}")

data = []
for line in lines:
    parts = line.strip().split()
    values = [float(x) for x in parts]
    data.append(values)

print(f"Data rows: {len(data)}")
if len(data) > 0:
    print(f"First row length: {len(data[0])}")
    print(f"First row: {data[0]}")
