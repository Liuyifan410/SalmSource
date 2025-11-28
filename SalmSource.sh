#!/usr/bin/env bash
#
#echo original parameters=[$@]
#
ARGS=`getopt -o hp:i:s: --long prefix:,input:,stage: -n "$0" -- "$@"`
if [ $? != 0 ]; then
    echo "Terminating..."
    exit 1
fi
#
#echo ARGS=[$ARGS]
#将规范化后的命令行参数分配至位置参数（$1,$2,...)
eval set -- "${ARGS}"
#echo formatted parameters=[$@]

while true
do
    case "$1" in
	-h|--help)
	    echo "-p|--prefix the prefix of the output file"
	    echo "-i|--input the input genome file, the input genome file should be fasta format"
	    echo "-s|--stage the stage of model, two options can be chosen: one or two"
      exit 0
	    ;;
        -p|--prefix)
	    forprefix=$2
            shift 2
            ;;
        -i|--input)
	    forinput=$2
            shift 2
            ;;
        -s|--stage)
	    forstage=$2
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Internal error!"
            exit 1
            ;;
    esac
done
#
echo $forprefix,$forinput,$forstage
#
#Rscript CorePanCompare.R 666_sample /home/micro/LZP/PanCoreSample/VCsample/VCtree/parsnp.tree
#Rscript ks_testff.R $forprefix $fortreefile $forpanfile $forcutoff
#
if [ "$forstage" = "one" ]; then
	Rscript model1.R ${forinput} ${forprefix}
elif [ "$forstage" = "two" ]; then
	Rscript model2.R ${forinput} ${forprefix}
else
    echo "error, you should choose one or two as input!" >&2
    exit 1  # 终止脚本并返回错误状态
fi
