Thanks for the report.

This error indicates that the SQLite file passed to `MeSHDbi::MeSHDb()` does not contain the `DATA` table expected by an organism-specific MeSHDb annotation resource:

```text
no such table: DATA
```

The likely cause is the example line:

```r
filepath_hsa <- qr_hsa[[1]]
```

`AnnotationHub::query()` can return multiple resources, and the first result is not guaranteed to be the organism-specific MeSHDb annotation database required by `meshdata()`. The query result should be inspected first, then the appropriate Homo sapiens MeSHDb resource should be selected.

I updated the package code to validate `MeSHDb` objects before querying them. With the update, incompatible AnnotationHub resources will produce a clearer error explaining that the object must contain the required `DATA` columns (`GENEID`, `MESHCATEGORY`, `MESHID`, `SOURCEDB`) instead of failing with the low-level SQLite message.
