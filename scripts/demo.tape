# scripts/demo.tape
# VHS script to record cellify usage

Output docs/images/demo.gif

Set FontSize 14
Set FontFamily "Monaco"
Set Width 800
Set Height 520
Set WindowBar Colorful
Set Theme "Dracula"
Set Shell zsh
Set TypingSpeed 60ms

# Start recording
Sleep 1s

# Clear terminal screen
Type "clear" Enter
Sleep 500ms

# Activate virtual environment
Type "source .venv/bin/activate" Enter
Sleep 1s

# Query TiO2
Type "cellify -i TiO2" Enter
Sleep 3.5s

# Select structure candidate 1 (COD 1544349)
Type "1" Enter
Sleep 4.5s

# Run conventional cell conversion and 2x2x4 supercell generation
Type "cellify -i TiO2_1544349_supercell.cif --conventional -d 2 2 4 -o TiO2_supercell.cif" Enter
Sleep 5s

# Clean up generated files to keep the workspace clean
Type "rm TiO2_1544349_supercell.cif TiO2_supercell.cif" Enter
Sleep 1s

# Final pause
Sleep 2s
