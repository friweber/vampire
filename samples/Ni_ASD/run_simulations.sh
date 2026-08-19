#!/bin/bash
#------------------------------------------------------------------
# Run Ni ASD quantum Heun simulations for three noise types:
#   classical        — white Gaussian noise (standard LLG baseline)
#   quantum          — coth spectral weight with zero-point energy
#   quantum-no-zero  — (coth-1) thermal fluctuations only
#
# Usage: ./run_simulations.sh [path-to-vampire-serial]
# Default binary: ../../vampire-serial
#------------------------------------------------------------------

VAMPIRE=${1:-../../vampire-serial}

if [ ! -f "$VAMPIRE" ]; then
    echo "Error: vampire binary not found at $VAMPIRE"
    echo "Usage: $0 [path-to-vampire-serial]"
    exit 1
fi

CASES=( "classical" "quantum" "quantum_no_zero" )
INPUTS=( "input_classical" "input_quantum" "input_quantum_no_zero" )

for i in "${!CASES[@]}"; do
    CASE="${CASES[$i]}"
    INPUT="${INPUTS[$i]}"
    OUTDIR="results_${CASE}"

    echo "========================================"
    echo "Running: $CASE"
    echo "Input:   $INPUT"
    echo "Output:  $OUTDIR/"
    echo "========================================"

    mkdir -p "$OUTDIR"

    # Run vampire from output directory so all output files land there.
    # Symlink the mat file so vampire can find it.
    cp "$INPUT" "$OUTDIR/input"
    cp "Ni_ASD.mat" "$OUTDIR/Ni_ASD.mat"

    (cd "$OUTDIR" && "../$VAMPIRE" > vampire.log 2>&1)

    if [ $? -eq 0 ]; then
        echo "Completed successfully. Output in $OUTDIR/"
    else
        echo "ERROR: simulation failed — check $OUTDIR/vampire.log"
    fi
    echo ""
done

echo "All simulations complete."
echo "Output directories: results_classical/  results_quantum/  results_quantum_no_zero/"
