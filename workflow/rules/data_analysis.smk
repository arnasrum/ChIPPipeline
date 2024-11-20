

rule compute_matrix:
    #input:
        #"test"
    output:
        "test"
    conda:
        "../envs/data_analysis.yml"
    params:
        ""
    shell:
        """
        touch test 
        """
