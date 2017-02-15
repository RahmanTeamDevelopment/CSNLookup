# CSNLookup

A web application written in Flask available online in web browsers (so it is easily accessible in the lab).

The user can submit a list of annotations and the webpage will tell if the submitted annotations are of correct CSN v1.0 format, and whether they correspond to potential variants in the provided transcripts. If there is any issue with an annotation, the tool reports a detailed explanation why the annotation has failed. The tool can detect several types of problems from issues in the syntax of CSN to more complex things such as incorrect indel alignment etc. 

If the annotation is correct, details of the variant (e.g. genome coordinates, alleles, location, class) are provided.

There is also an interactive page "Transcript catalog" and a "User manual" page included.

The tool will save all submissions (input and output files) on the server for possible later interest.

See the file INSTALL for steps of installation before deployment to server.



