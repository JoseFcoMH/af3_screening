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
        f"{ENV_DIR}/base.yaml"
    params:
        prediction_dirs=lambda wc, input: [
            str(Path(model).parent) for model in input.models
        ],
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail

        python {SCRIPT_DIR}/project_to_structure.py \
            --af3_list {params.prediction_dirs:q} \
            --out_file {output.attr_file:q} \
            --cif_file {input.cif_structure:q} \
            --n_cpus {threads} \
            --reference_chain {wildcards.prey}
        """


rule apply_defattr:
    input:
        cif=lambda wc: BAITS_DICT[wc.bait],
        defattr="<pair_output>/{bait}__{prey}.defattr",
    output:
        cxs="<pair_output>/{bait}__{prey}.cxs",
    log:
        "<logs>/apply_defattr/<pair_name>.log",
    params:
        modules=modules_for("apply_defattr"),
        attribute_name=lambda wc: f"{wc.bait}_{wc.prey}",
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail

        {params.modules}

        CHIMERAX=$(command -v chimerax || command -v ChimeraX)

        "$CHIMERAX" --nogui <<'EOF'
open {input.cif}
open {input.defattr}
color byattribute r:iCP_{params.attribute_name} #!1 target csab palette 0,#0e8bff:0.5,white:1,red noValueColor #999999
save {output.cxs}
exit
EOF
        """
