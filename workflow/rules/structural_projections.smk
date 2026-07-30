import os

UDOCKER_DIR = Path(
    config.get("udocker_dir", os.environ.get("UDOCKER_DIR", "~/.udocker"))
)
CHIMERAX_IMAGE = config.get("chimerax_image", "maduprey/chimerax:1.5")
CHIMERAX_NAME = config.get(
    "chimerax_container_name", CHIMERAX_IMAGE.split("/")[-1].replace(":", "-")
)


rule setup_chimerax:
    output:
        directory(UDOCKER_DIR / "containers" / f"{CHIMERAX_NAME}.ready"),
    log:
        "<logs>/setup_chimerax.log",
    conda:
        ENVS["udocker"]
    params:
        udocker_dir=str(UDOCKER_DIR),
        chimerax_image=CHIMERAX_IMAGE,
        chimerax_container=CHIMERAX_NAME,
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail

        export UDOCKER_DIR={params.udocker_dir:q}

        mkdir -p "$UDOCKER_DIR"

        udocker install

        if ! udocker ps | grep -Fxq {params.chimerax_container:q}; then
            udocker pull {params.chimerax_image:q}
            udocker create \
                --name={params.chimerax_container:q} \
                {params.chimerax_image:q}
            udocker setup --execmode=P1 {params.chimerax_container:q}
        fi

        touch {output:q}
        """


rule map_prediction_to_structure:
    input:
        confs=confidences_for_pair,
        models=models_for_pair,
        cif_structure=lambda wc: BAITS_DICT[wc.bait],
    output:
        attr_file="<pair_output>/{bait}__{prey}.defattr",
    log:
        "<logs>/map_prediction_to_structure/<pair_name>.log",
    conda:
        ENVS["base"]
    params:
        prediction_dirs=lambda wc, input: [
            str(Path(model).parent) for model in input.models
        ],
        script_dir=SCRIPT_DIR,
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail

        python {params.script_dir}/project_to_structure.py \
            --af3_list {params.prediction_dirs:q} \
            --out_file {output.attr_file:q} \
            --cif_file {input.cif_structure:q} \
            --n_cpus {threads} \
            --reference_chain {wildcards.prey}
        """


rule apply_defattr:
    input:
        container=ancient(rules.setup_chimerax.output),
        cif=lambda wc: BAITS_DICT[wc.bait],
        defattr="<pair_output>/{bait}__{prey}.defattr",
    output:
        cxs="<pair_output>/{bait}__{prey}.cxs",
    log:
        "<logs>/apply_defattr/<pair_name>.log",
    conda:
        ENVS["udocker"]
    params:
        attribute_name=lambda wc: f"{wc.bait}_{wc.prey}",
        cif_dir=lambda wc, input: str(Path(input.cif).parent),
        defattr_dir=lambda wc, input: str(Path(input.defattr).parent),
        udocker_dir=UDOCKER_DIR,
        udocker_container=CHIMERAX_NAME,
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail

        export UDOCKER_DIR={params.udocker_dir:q}

        udocker run \
            --nobanner \
            --volume="$PWD:$PWD" \
            --volume={params.cif_dir}:{params.cif_dir} \
            --volume={params.defattr_dir}:{params.defattr_dir} \
            --workdir="$PWD" \
            {params.udocker_container:q} \
            ChimeraX --nogui <<'EOF'
open {input.cif}
open {input.defattr}
color byattribute r:iCP_{params.attribute_name} #!1 target csab palette 0,#0e8bff:0.5,white:1,red noValueColor #999999
save {output.cxs}
exit
EOF
        """
