# Error Contract 1.0

Error codes are stable machine interfaces. Diagnostic messages may gain detail
and MUST NOT be parsed. A failed analysis returns no partial scientific JSON.

## Native analysis codes

| Code | Meaning | Retry guidance |
| --- | --- | --- |
| `MITO-E1001` | invalid configuration or API argument | correct configuration |
| `MITO-E1101` | input cannot be opened | correct path/permissions |
| `MITO-E1102` | unsupported input format | convert/use a supported format |
| `MITO-E1103` | malformed input record | repair or regenerate input |
| `MITO-E1104` | input contains no records | provide non-empty input |
| `MITO-E1201` | reference cannot be opened | correct path/permissions |
| `MITO-E1202` | invalid reference content | provide a valid single FASTA record |
| `MITO-E1301` | required resource cannot be opened | restore the pinned bundle |
| `MITO-E1302` | required resource is malformed/incomplete | replace with verified bundle |
| `MITO-E1401` | compiled dependency/capability unavailable | install dependency and rebuild |
| `MITO-E1501` | cooperative cancellation | submit a new job if desired |
| `MITO-E1601` | memory allocation or configured observation/phase/resource bound was exceeded | reduce/page the workload or increase a reviewed limit |
| `MITO-E9001` | internal invariant or unknown failure | preserve evidence and report defect |

The C ABI exposes `mito_engine_get_last_error_code()` and
`mito_engine_get_last_error()` as thread-local state, plus
`mito_engine_error_schema_version()`. Rust returns `mito_ffi::MitoError` with
public `code` and `message` fields.

## HTTP error envelope

Every synchronous API error uses:

```json
{
  "error": {
    "schema_version": "1.0",
    "code": "MITO-API-E1005",
    "message": "upload exceeds the configured byte limit",
    "retryable": false
  }
}
```

HTTP-level codes cover invalid request (`E1001`), authorization (`E1002`), not
found (`E1003`), conflict (`E1004`), payload limit (`E1005`), capacity (`E1006`),
not ready (`E1007`), service unavailable (`E2001`), and internal failure
(`E9001`). Completed job failures preserve the native `MITO-E...` code in the
same envelope. Internal filesystem/process details are logged server-side and
are not returned to clients.

Adding a code is backward compatible. Reassigning a code or changing its
meaning requires a new error-schema version.
