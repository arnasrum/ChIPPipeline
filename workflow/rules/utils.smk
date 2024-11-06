
rule unzip:
    input:
        "{file}.gz"
    output:
        "{file}"
    wildcard_constraints:
        file = r"^([\/A-Za-z0-9])*"
    shell:
        "gzip -dk {input} -f"