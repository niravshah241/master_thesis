#!/bin/bash

RELDIR=`dirname $0`

[ -f $RELDIR/mtocpp.conf ] && cd $RELDIR

if [ ! -f ./mtocpp.conf ]; then
    echo "Error: make_docu.sh needs to be run in project directory!"
    exit 1;
fi

if [ ! -f ./doxygen/configuration ]; then
    echo "configuration file does not exist. Shall it be constructed? (y/n)"
    while [ "z$answer" != "zn" -a "z$answer" != "zy" ]; do
        read answer
    done
    [ $answer == "n" ] && exit 2;

    echo
    mkdir -p ./doxygen/

    BASEDIR_GUESS=${PWD}
    echo "Please give me the base directory of the RBmatlab installation."
    echo "empty line for default = [ $BASEDIR_GUESS ]"
    read BASEDIR
    [ -z "$BASEDIR" ] && BASEDIR=$BASEDIR_GUESS

    echo 'define(`BASEDIR'\'', `'"$BASEDIR"\'')' > ./doxygen/configuration

    # Try to find a working mtocpp installation:
    MTOCPP_PATH="$(which mtocpp)"
    if [ -n "${MTOCPP_PATH}" ]; then
        MTOCPP_DIR=$(dirname $MTOCPP_PATH)
        echo "Using mtocpp installed in ${MTOCPP_DIR}!"
    else
        if [ -d "$BASEDIR/3rdparty/mtocpp" ];
        then
            MTOCPP_DIR=$BASEDIR/3rdparty/mtocpp
        elif [ -f "$BASEDIR/3rdparty/mtocpp" ];
        then
            MTOCPP_DIR=$BASEDIR/3rdparty
        else
            echo "Could not find mtoc directory in \"$BASEDIR/3rdparty/mtocpp\""
            echo "Please give a valid destination for mtocpp filter:"
            read MTOCPP_DIR
        fi
    fi

    MTOCPP_DIR=$MTOCPP_DIR/

    MTOCPP_SUFFIX=""
    SUFFIX_LIST="i386 i686 x86_64 win32 win64 mac64 none"
    for MTOCPP_SUFFIX in " " $SUFFIX_LIST
    do
        CANDIDATE=${MTOCPP_DIR}mtocpp$MTOCPP_SUFFIX
        echo $CANDIDATE
        [ -x $CANDIDATE ] && FIRST_LINE=$($CANDIDATE --help | head -1)
        if echo $FIRST_LINE | grep -q "Usage"
        then
            echo "Found old version of mtocpp"
            break
        elif echo $FIRST_LINE | grep -q "Version"
        then
            break
        fi
    done

    if [ "z$MTOCPP_SUFFIX" == "znone" ]
    then
        echo "Error: Could not find a working mtocpp filter."
        echo "Please install it first for example by downloading from"
        echo "http://www.morepas.org/software/mtocpp/"
        exit;
    else
        echo "Found mtocpp $CANDIDATE and it works".
    fi


    echo 'define(`MTOCPP_DIR'\'', `'"$MTOCPP_DIR"\'')' >> ./doxygen/configuration
    echo 'define(`MTOCPP_SUFFIX'\'', `'"$MTOCPP_SUFFIX"\'')' >> ./doxygen/configuration
fi

LATEX=0
if [ $# -eq 1 ] && [ $1 == "latex" ];
then
    sed -i -e 's/LATEX_OUTPUT := [^;]*;/LATEX_OUTPUT := true;/' mtocpp.conf
    if ! grep -q LATEX_OUTPUT mtocpp.conf;
    then
        sed -i -e '1 s/^/LATEX_OUTPUT := true;\n/' mtocpp.conf
    fi
    m4 ./doxygen/Doxyfile_latex.in > ./doxygen/Doxyfile_latex
    LATEX=1;
else
    sed -i -e 's/LATEX_OUTPUT := [^;]*;/LATEX_OUTPUT := false;/' mtocpp.conf
    m4 ./doxygen/Doxyfile.in > ./doxygen/Doxyfile

fi

TOEVAL=$(echo -e "include(\`./doxygen/configuration')SUFFIX=\"MTOCPP_SUFFIX\"\nEXEC_DIR=\"MTOCPP_DIR\"\n" | m4 -)
eval "${TOEVAL}"


if [ $LATEX -ne 0 ]; then
    doxygen doxygen/Doxyfile_latex 2>&1 1> doxygen/doxygen_latex.log | grep -v synupdate | grep -v docupdate | grep -v display | grep -v subsref | tee doxygen/doxygen_latex.err
    cd doxygen/latex

    for i in *.tex; do ${MTOCPP_DIR}mtocpp_post${SUFFIX} $i; done
else
    doxygen doxygen/Doxyfile 2>&1 1>doxygen/doxygen.log \
        | tee doxygen/doxygen.err | grep -v synupdate \
        | grep -v docupdate \
        | grep -v display \
        | grep -v subsref \
        | grep -v "is not documented" \
        | tee doxygen/doxygen_small.err

    grep "is not documented" doxygen/doxygen.err \
        | perl -ne '@_=split(/:/); print "$_[0]\n"' \
        | sort \
        | uniq > doxygen/needs_docu.log

    N_UNDOC="$(wc -l doxygen/needs_docu.log | cut -f 1 -d ' ')"
    N_ERRORS="$(wc -l doxygen/doxygen.err | cut -f 1 -d ' ')"
    echo "$N_UNDOC files have incomplete doxygen documentation."
    echo "See doxygen/doxygen.err and doxygen/needs_docu.log for more details."
    echo
    echo "$N_ERRORS warnings or errors are thrown by doxygen."
    echo "See doxygen/doxygen.err for more details."

    ${EXEC_DIR}mtocpp_post${SUFFIX} doxygen/html;

fi

