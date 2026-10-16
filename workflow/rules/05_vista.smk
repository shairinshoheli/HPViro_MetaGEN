# ============================================================
# Step 5: VISTA metagenomic subspecies (mgSs) + community state types (mgCST)
# Input: the compiled VIRGO2 summary.NR table (step 4).
# ============================================================

rule install_vista:
    output:
        touch(VISTA_READY),
    log:
        f"{LOG_DIR}/install_vista.log",
    shell:
        r"""
        set -euo pipefail
        exec 3>&2                     # keep a handle on the terminal for error messages
        exec > {log} 2>&1

        rm -rf "{VISTA_HOME}"
        git clone https://github.com/JHolm-Lab/VISTA.git "{VISTA_HOME}"

        # Use a manually downloaded archive if present, otherwise try to download it
        if [ ! -s "{VISTA_ARCHIVE}" ]; then
            rm -f "{VISTA_ARCHIVE}"
            curl -fL --retry 3 -o "{VISTA_ARCHIVE}" "{VISTA_DATA_URL}" || true
        fi

        if [ ! -s "{VISTA_ARCHIVE}" ] || ! tar -tzf "{VISTA_ARCHIVE}" > /dev/null 2>&1; then
            rm -f "{VISTA_ARCHIVE}"
            echo "ERROR: VISTA data archive is missing or not a valid .tar.gz." >&3
            echo "Download VISTA_data.tar.gz in a browser from https://doi.org/10.6084/m9.figshare.28684934.v10" >&3
            echo "save it as {VISTA_ARCHIVE}, then run again (details: {log})." >&3
            exit 1
        fi

        tar -xzf "{VISTA_ARCHIVE}" -C "{VISTA_HOME}"

        # Fail loudly if the classifier script or the data folder is missing
        test -f "{VISTA_HOME}/run_VISTA.R"
        test -d "{VISTA_HOME}/VISTA_data"
        """


# run_VISTA.R writes timestamped files to the current directory, so we cd in.
rule vista:
    input:
        nr=f"{VIRGO2_COMP_DIR}/virgo2_compiled.summary.NR.txt",
        ready=VISTA_READY,
    output:
        directory(VISTA_DIR),
    log:
        f"{LOG_DIR}/vista.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{output}"
        in=$(realpath {input.nr})
        vista=$(realpath {VISTA_HOME})
        logfile=$(realpath -m {log})

        cd "{output}"
        Rscript "$vista/run_VISTA.R" "$in" "$vista" > "$logfile" 2>&1

        # Fail loudly if the classification table was not written
        compgen -G "mgCSTs_*.csv" > /dev/null
        """
