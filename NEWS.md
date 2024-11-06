# lfe 3.1.0

* Added a `NEWS.md` file to track changes to the package.
* The use of NAMED and SET_NAMED in R C code is not recommended as they are not
  part of the public API and their behavior is subject to change. This led to
  change `utils.c`.
* Uses DOI numbers instead of URLs for references in the documentation.
