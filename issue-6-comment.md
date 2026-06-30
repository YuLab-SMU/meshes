Thanks for the report.

The original `trying to get slot "result" from an object of a basic class ("NULL")` error happens after the enrichment engine returns `NULL` because none of the input Entrez Gene IDs can be mapped to the selected MeSH annotation set.

This is usually caused by one of the following:

- the selected `MeSHDb` resource is not the expected organism-specific annotation database;
- the selected `database` / `category` combination has no overlap with the input genes;
- the example input gene IDs are not represented in the selected MeSH annotation source.

I checked the related code path while updating the package. The current code no longer assumes the enrichment result is always non-`NULL`, and I also added validation for the supplied `MeSHDb` resource so that incompatible AnnotationHub results fail earlier with a clearer message. For examples and tutorials, we will keep the runnable workflows in the online book rather than relying on a brittle `qr_hsa[[1]]` selection from AnnotationHub.
