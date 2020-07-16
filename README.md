# Exciting Parser

This is the parser for [exciting](http://exciting-code.org/).
It is part of the [NOMAD Laboratory](http://nomad-lab.eu).
The official version lives at

    git@gitlab.mpcdf.mpg.de:nomad-lab/parser-exciting.git

you can browse it at

    https://gitlab.mpcdf.mpg.de/nomad-lab/parser-exciting

It relies on having the nomad-meta-info and the python common repositories one level higher.
The simplest way to have this is to check out nomad-lab-base recursively:

    git clone --recursive git@gitlab.mpcdf.mpg.de:nomad-lab/nomad-lab-base.git

then this will be in parsers/exciting.


## Input Filenames

Variable names & variable definitions are case sensitive

Some variable names point to other variable names, e.g., `inputgwFile`. <br>

Sorted by variable name

ToDo:
- split the table into two groups: necessary and optional input files.
- search for changes of `PWD`:  `grep -H -r  "chdir"` # os.chdir()
- A "developers" version of the filename list could look like the table below, however, the "users" version should just be a concise version of the second column, using examples instead of nested definitions. E.g.,using `EPSILON_*.OUT` instead of `'EPSILON_' + ext + ... +.'OUT'`

|VARIABLE NAME |      DEFINITION|
|--- | --- |
|'INFO_VOL'       |    --- |
|'str.out'        |    ---|
|DielNoSymFile    | 'DIELTENS0_NOSYM' + qExt00 + '.OUT' |
|DielSymFile     | 'DIELTENS0' + qExt00 + '.OUT'|
|QFile           | "QPOINTS.OUT"|
|bandBorGWFile   | "BAND-QP.OUT" |
|bandCarbGWFile | "bandstructure-qp.dat" |
|bandFile |  "bandstructure.xml" |
|dosFile | "dos.xml"|
|dosGWFile | "TDOS-QP.OUT"|
|eigvalFile | "EIGVAL.OUT"|
|eigvalGWFile | "EVALQP.DAT" or "EVALQP.TXT"|
|**epsFile** | ???? |
|epsilonLocalField  | 'EPSILON_' + ext + 'FXC' + self.tddftKernel[0] + '_OC' + tensorComp[j] + qExt00 + '.OUT' |
|epsilonNoLocalField  | 'EPSILON_' + ext + 'NLF_' + 'FXC' + self.tddftKernel[0] + '_OC' + tensorComp[j] + qExt00 + '.OUT'|
|**excFile** | ?????? |
|fermiSurfFile  | "FERMISURF.bxsf"|
|gFile -> gw_file | "GW_INFO.OUT"|
|inputGSFile  |"input.xml" |
|inputXSFile  | "input.xml"|
|inputgwFile = [inputgw1File, inputgw2File, inputFile] | ["input-gw.xml" , "input.xml"] |
|lossFunctionLocalFieldFile  | 'LOSS_' + ext + 'FXC' + self.tddftKernel[0] + '_OC' + tensorComp[j] + qExt00 + '.OUT'|
|lossFunctionNoLocalFieldFile  |'LOSS_' + ext + 'NLF_' + 'FXC' + self.tddftKernel[0] + '_OC' + tensorComp[j] + qExt00 + '.OUT'|
|outputEpsFile  |"EPSILON_BSE" + self.bsetype + '_SCR' + self.screentype + "_OC" + self.tensorComp[i] + ".OUT" |
|outputSigmaFile  |"SIGMA_BSE" + self.bsetype + '_SCR' + self.screentype + "_OC" + self.tensorComp[i] + ".OUT" |
|outputXSFile  | "EXCITON_BSE" + self.bsetype + '_SCR' + self.screentype + "_OC" + self.tensorComp[i] + ".OUT"|
|qPlusGFile  |'GQPOINTS' + qExt00 + '.OUT'|
|sigmaLocalFieldFile  | 'SIGMA_' + ext + 'FXC' + self.tddftKernel[0] + '_OC' + tensorComp[j] + qExt00 + '.OUT'|
|sigmaNoLocalFieldFile |'SIGMA_' + ext + 'NLF_' + 'FXC' + self.tddftKernel[0] + '_OC' + tensorComp[j] + qExt00 + '.OUT' |
|vertexGWFile  | "BANDLINES.OUT"|
|vertexLabGWFile  | "bandstructure.xml"|

