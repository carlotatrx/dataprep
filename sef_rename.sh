#!/bin/bash

# --- CONFIGURACIÓ ---
# Ruta absoluta on hi ha les dades
TARGET_DIR="/scratch3/PALAEO-RA/daily_data/final"
# Fitxer on guardarem el registre (es crearà on executis l'script)
LOG_FILE="output.txt"

echo "Inici del procés - $(date)" > "$LOG_FILE"
echo "Directori objectiu: $TARGET_DIR" >> "$LOG_FILE"
echo "------------------------------------------" >> "$LOG_FILE"

# Comprovem que el directori existeix abans de començar
if [ ! -d "$TARGET_DIR" ]; then
    echo "[ERROR CRÍTIC] El directori $TARGET_DIR no existeix." | tee -a "$LOG_FILE"
    exit 1
fi

# Busquem recursivament fitxers .tsv DINS del directori objectiu
find "$TARGET_DIR" -type f -name "*.tsv" -print0 | while IFS= read -r -d '' file; do
    
    # Obtenim noms i rutes
    dir_path=$(dirname "$file")
    folder_name=$(basename "$dir_path") # Nom de la carpeta (ex: Geneva)
    file_name=$(basename "$file")       # Nom del fitxer (ex: Geneva_1782...)

    # ---------------------------------------------------------
    # CONDICIÓ: El fitxer comença pel nom de la carpeta?
    # ---------------------------------------------------------
    if [[ "$file_name" != "$folder_name"* ]]; then
        # Si NO comença pel nom de la carpeta, assumim que ja té prefix.
        echo "[SKIP] Ja té nom de projecte: $file_name (Carpeta: $folder_name)" >> "$LOG_FILE"
        continue
    fi

    # ---------------------------------------------------------
    # EXTRACCIÓ: Busquem el 'Source'
    # ---------------------------------------------------------
    # Llegim la línia 'Source' i traiem el retorn de carro (\r) per si ve de Windows
    project_name=$(awk -F'\t' '$1 == "Source" {print $2}' "$file" | tr -d '\r')

    # Validació: Si 'Source' està buit
    if [ -z "$project_name" ]; then
        echo "[ERROR] No s'ha trobat 'Source' a: $file" >> "$LOG_FILE"
        continue
    fi

    # Netegem espais en blanc sobrers
    project_name=$(echo "$project_name" | xargs)

    # ---------------------------------------------------------
    # CANVI DE NOM
    # ---------------------------------------------------------
    new_base_name="${project_name}_${file_name}"
    new_full_path="${dir_path}/${new_base_name}"

    mv "$file" "$new_full_path"

    if [ $? -eq 0 ]; then
        echo "[OK] Renomenat: $file_name  -->  $new_base_name" >> "$LOG_FILE"
    else
        echo "[FAIL] Error movent: $file" >> "$LOG_FILE"
    fi

done

echo "------------------------------------------" >> "$LOG_FILE"
echo "Procés finalitzat." >> "$LOG_FILE"

echo "Feta la feina. Revisa $LOG_FILE per veure els resultats."