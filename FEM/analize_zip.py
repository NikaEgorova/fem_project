import csv
import os
import re
import zipfile

# Назва вашого архіву
ZIP_FILE_PATH = "FEM_fortran.zip"
OUTPUT_CSV = "fortran_analysis.csv"


def clean_fortran_line(line):
    """Видаляє коментарі з рядків Fortran (Fixed & Free format)"""
    # Коментарі у Fixed format (літери C, c, * у першій колонці)
    if len(line) > 0 and line[0] in ["C", "c", "*"]:
        return ""
    # Коментарі знак '!'
    if "!" in line:
        line = line.split("!")[0]
    return line.strip()


def analyze_fortran_code(zip_path, output_csv):
    results = []

    # Регулярні вирази для пошуку оголошень та викликів
    def_pattern = re.compile(
        r"^\s*(?:[\w\*]+\s+)*(?:subroutine|function)\s+([a-zA-Z0-9_]+)",
        re.IGNORECASE,
    )
    call_pattern = re.compile(r"\bcall\s+([a-zA-Z0-9_]+)", re.IGNORECASE)

    with zipfile.ZipFile(zip_path, "r") as z:
        # Фільтруємо лише файли вихідного коду (.for, .f, .FOR, .F)
        fortran_files = [
            f
            for f in z.namelist()
            if f.lower().endswith((".for", ".f")) and not f.startswith(".git/")
        ]

        for file_name in fortran_files:
            try:
                with z.open(file_name) as f:
                    content = f.read().decode("utf-8", errors="ignore")
            except Exception:
                continue

            lines = content.splitlines()

            defined_subprograms = set()
            called_subprograms = set()

            for line in lines:
                cleaned = clean_fortran_line(line)
                if not cleaned:
                    continue

                # Пошук оголошених підпрограм у цьому файлі
                def_match = def_pattern.match(cleaned)
                if def_match:
                    defined_subprograms.add(def_match.group(1).lower())

                # Пошук викликів підпрограм (CALL name)
                calls = call_pattern.findall(cleaned)
                for c in calls:
                    called_subprograms.add(c.lower())

            # Розподіл викликів на внутрішні та зовнішні
            internal_calls = called_subprograms.intersection(
                defined_subprograms
            )
            external_calls = called_subprograms - defined_subprograms

            # Очищення шляху файлу для зручності
            short_name = os.path.basename(file_name)

            results.append({
                "filename": short_name,
                "external_calls": ", ".join(sorted(external_calls))
                if external_calls
                else "-",
                "internal_calls": ", ".join(sorted(internal_calls))
                if internal_calls
                else "-",
            })

    # Запис у CSV файл
    with open(output_csv, mode="w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f, delimiter=";")
        # Заголовки стовпчиків
        writer.writerow([
            "Назва файлу",
            "Назви сторонніх підпрограм",
            "Назви внутрішніх підпрограм",
        ])

        for row in results:
            writer.writerow([
                row["filename"],
                row["external_calls"],
                row["internal_calls"],
            ])

    print(
        f"Аналіз завершено! Файл '{output_csv}' успішно створено для Google Calc."
    )


if __name__ == "__main__":
    analyze_fortran_code(ZIP_FILE_PATH, OUTPUT_CSV)