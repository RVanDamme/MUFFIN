# Seqera MCP connection and discovery

Use the plugin's configured Seqera MCP connection (`https://mcp.seqera.io/mcp`).
The host manages OAuth sign-in and credentials. Never paste tokens into prompts,
files or command examples. A missing connection is a capability gap, not a reason
to bypass authentication with curl or copied credentials.

## Discover → execute → verify

1. Inspect the tools actually exposed by the host. Their current descriptions
   own API names, parameter schemas and authentication behavior.
2. Use `search_seqera_api` with the operation's intent. Read the returned full
   schema, including optional fields, nesting and types; suggestions alone may
   omit optional launch settings. Resolve workspace/resource IDs from returned
   metadata rather than guessing.
3. Use `call_seqera_api` with the discovered service, exact `api_name` and
   schema-shaped parameters. This is not a raw HTTP method/endpoint/body API.
4. For mutations, confirm the authorized scope and verify the resulting state
   with an authoritative read. Broaden discovery or report a capability gap when
   an operation is unavailable; do not fabricate names or mirror API schemas.

For resume, discover the existing-workflow resume operation and preserve the
original launch/session relationship. A Platform workflow ID is not a Nextflow
session UUID. For saved-pipeline or ad-hoc launch, follow `launch-workflow`'s
relevant procedure and preserve the selected revision, parameters and lineage.

Use exposed nf-core discovery helpers for module selection and inspect their
input/output contracts. Workspace execution, resource access and artifact
handoff remain host responsibilities; a tool response is not proof that a
container or pipeline ran successfully.
