parola=$1
for file in *; do
    # Verifica se il file è un file regolare e non una directory
    if [ -f "$file" ]; then
        # Controlla se la parola è presente nel file
        if grep -q "$parola" "$file"; then
            echo "La parola '$parola' è presente nel file: $file"
            grep --color -n "$parola" "$file"
        fi
    fi
done