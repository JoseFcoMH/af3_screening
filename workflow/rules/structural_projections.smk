import os

UDOCKER_DIR = Path(
    config.get("udocker_dir", os.environ.get("UDOCKER_DIR", "~/.udocker"))
)
CHIMERAX_IMAGE = config.get("chimerax_image", "maduprey/chimerax:1.5")
CHIMERAX_NAME = config.get("chimerax_container_name", "chimerax")


rule setup_chimerax:
    output:
        str(UDOCKER_DIR / "containers" / f"{CHIMERAX_NAME}.ready"),
    log:
        "<logs>/setup_chimerax.log",
    conda:
        ENVS["udocker"]
    params:
        udocker_dir=str(UDOCKER_DIR),
        image=CHIMERAX_IMAGE,
        container=CHIMERAX_NAME,
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail

        export UDOCKER_DIR={params.udocker_dir:q}

        mkdir -p "$UDOCKER_DIR"
        mkdir -p "$(dirname {output:q})"

        udocker install

        if ! udocker inspect {params.container:q} >/dev/null 2>&1; then
            udocker pull {params.image:q}

            udocker create \
                --name={params.container:q} \
                {params.image:q}
        fi

        udocker setup \
            --execmode=P1 \
            {params.container:q}

        # Do not create the sentinel unless the container really runs.
        udocker run \
            --nobanner \
            {params.container:q} \
            /bin/sh -c 'echo "ChimeraX container works"'

        touch {output:q}
        """


rule map_prediction_to_structure:
    input:
        confs=rules.get_pair_confidences.output.out_dir,
        models=rules.get_pair_models.output.out_dir,
        cif_structure=lambda wc: BAITS_DICT[wc.bait],
    output:
        attr_file="<pair_output>/{bait}__{prey}.defattr",
    log:
        "<logs>/map_prediction_to_structure/<pair_name>.log",
    conda:
        ENVS["base"]
    params:
        predictions=lambda wc, input: [
            Path(model).stem.split("_model.cif")[0] for model in input.models
        ],
        mock_dir=lambda wc, input: f"{SCRATCH}/.rundata/{wc.bait}__{wc.prey}/map_prediction_to_structure",
        mock_pred_dirs=lambda wc, input: [
            f"{SCRATCH}/.rundata/{wc.bait}__{wc.prey}/map_prediction_to_structure/{Path(model).stem.split('_model.cif')[0]}"
            for model in input.models
        ],
        script_dir=SCRIPT_DIR,
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail

        for pred in {params.predictions:q}; do
            mkdir -p "{params.mock_dir}/$pred"
            ln -s "{params.confs}/${{pred}}_model.cif" "{params.mock_dir}/${{pred}}/${{pred}}_model.cif"
        done

        python {params.script_dir}/project_to_structure.py \
            --af3_list {params.mock_pred_dirs:q} \
            --out_file {output.attr_file:q} \
            --cif_file {input.cif_structure:q} \
            --n_cpus {threads} \
            --reference_chain {wildcards.prey}

        rm -rf {params.mock_dir:q}
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
        chimerax_home=lambda wc: str(
            Path(SCRATCH) / f"chimerax_home/{wc.bait}__{wc.prey}"
        ),
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail

        export UDOCKER_DIR={params.udocker_dir:q}

        CHIMERAX_HOME={params.chimerax_home:q}
        mkdir -p "$CHIMERAX_HOME/.config/ChimeraX"

        udocker run \
            --nobanner \
            --nometa \
            --volume={params.cif_dir:q}:{params.cif_dir:q} \
            --volume={params.defattr_dir:q}:{params.defattr_dir:q} \
            --volume="$CHIMERAX_HOME:/home/chimerax" \
            --env=HOME=/home/chimerax \
            --env=XDG_CONFIG_HOME=/home/chimerax/.config \
            {params.udocker_container:q} \
            chimerax --nogui <<'EOF'
open {input.cif}
open {input.defattr}
color byattribute r:iCP_{params.attribute_name} #!1 target csab palette 0,#0e8bff:0.5,white:1,red noValueColor #999999
save {output.cxs}
exit
EOF
        """
