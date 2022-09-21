configfile: "./config.yaml"


rule build_bowtie_database:
    input:
        "code/build_database.sh"
    output: "logfile.log"
    threads: 4
    shell: 
       """ 
       bash {input}
       """