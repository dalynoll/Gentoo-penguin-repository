#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import os

# Directorios de entrada/salida
input_dir = "/data2/daly/gentoo/filogenias/orthofinder/1fasta32"  # CDS
output_dir = "/data2/daly/gentoo/filogenias/orthofinder/2fastaProt32"  # Proteínas
os.makedirs(output_dir, exist_ok=True)  # Crea el directorio si no existe

# Procesar cada archivo de CDS
for filename in os.listdir(input_dir):
    if filename.endswith(".fa") or filename.endswith(".fasta"):
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename.replace(".fa", "_prot.fa").replace(".fasta", "_prot.fasta"))

        print(f"Traduciendo {filename} → {output_path}")

        with open(output_path, "w") as output_handle:
            for record in SeqIO.parse(input_path, "fasta"):
                try:
                    protein_seq = Seq(record.seq).translate(to_stop=True)  # Traducción estándar sin STOPs
                    output_handle.write(f">{record.id}\n{protein_seq}\n")
                except Exception as e:
                    print(f"Error traduciendo {record.id}: {e}")

print("¡Traducción completada! Archivos guardados en:", output_dir)
