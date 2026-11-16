---
name: seqera-data-links
description: >
  Add, list, inspect, browse, update, and delete data links in Nextflow Platform
  workspaces. Use when asked to connect S3/GCS/Azure storage, manage Data Explorer
  entries, create data links, or troubleshoot bucket validation and storage access.
---
<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->

## Using this skill

Use the tools and runtimes available in the user's session. Establish the
pipeline, input data and execution environment from their request or current
context. Use user-selected locations for files and results. Consult companion
skills, bundled references and optional helpers by name when available.

For Seqera operations, read the [MCP connection](../launch-workflow/references/seqera-mcp/README.md) connection reference. Any `platform_*` examples
describe API operations. Discover their exact names and parameter schemas with
`search_seqera_api`, then invoke `call_seqera_api` using those schemas. The host
manages OAuth for the connected MCP server.


# Seqera Data Links

Manage Data Explorer entries in the user's existing Nextflow Platform workspace
through the production Seqera MCP connection. The client completes OAuth when
connecting; the MCP server handles Platform authentication. Never request,
print, or copy access tokens, bearer headers, cloud keys, or credential contents.

## Discover before calling

Read the [MCP connection](../launch-workflow/references/seqera-mcp/README.md) connection reference for discovery and authentication. Every
Platform operation follows:

1. Call `search_seqera_api` with the operation's intent, such as "list data links
   in a workspace", "create a data link", or "browse a data link".
2. Read the selected result's full parameter schema, including optional fields
   and nested objects. Required-field suggestions alone may omit needed options.
3. Call `call_seqera_api(service="platform", api_name=<returned api_name>,
   parameters=<arguments matching that schema>)`.

Names such as `platform_create_data_link` identify API operations discovered by
search; they are not separate tools to invoke directly. Use the exact API name,
field casing, nesting, and types returned by the server. Do not convert the
schema to snake_case or send raw HTTP payloads. Use MCP for bulk operations too.
If discovery cannot find an operation, broaden the search and report the
capability gap; do not bypass the connection with Platform curl or `tw` calls.

## Resolve the workspace and link

Use the workspace the user has already selected. If its ID is unknown, discover
workspace listing and retrieve accessible workspaces through `call_seqera_api`;
match the user's organization/workspace name to the returned ID. Ask for a
selection only when no workspace is specified or multiple matches remain.

Discover listing and lookup operations before changing an existing link. Use
returned link IDs and inspect the current name, provider, storage URI, description,
and status. Paginate listings when necessary; an incomplete page does not prove
a link is absent. Listing, lookup, and browsing are read-only operations.

## Create or update a link

For creation, establish:

- The target workspace and a descriptive link name.
- The exact user-selected storage URI and its provider.
- Whether the data is public or requires an existing workspace credential.
- Any description or other options requested by the user.

Discover the credential-listing operation when a private link needs credentials.
Return only identifying metadata needed to select an existing credential, such
as its ID, name, and provider; never expose credential secrets. If none matches,
use `ce-credentials-setup` for the separate credential setup workflow.

Discover the create operation and map these values to its current schema.
For example, a schema may call the URI `resourceRef`, require a bucket `type`,
or expect `credentialsId`; follow the returned schema rather than assuming those
fields. Private storage generally requires a compatible workspace credential.

For updates, discover the update operation and pass the returned link ID plus
only the requested changes. For deletion, resolve the exact link the user asked
to remove, then use the discovered delete operation. Do not delete or replace
other links while troubleshooting.

## Diagnose storage validation failures

Seqera validates the bucket or storage resource using the selected workspace
credential. A nonexistent bucket cannot be fixed by creating a Data Explorer
entry. A prefix inside a real bucket can be suitable even when that prefix has
not yet been populated; distinguish bucket existence from object availability.

For "bucket does not exist" or access errors, inspect the actual returned error
and the link's provider, URI, and credential metadata. Check for a typo, wrong
cloud account, wrong credential/provider, or insufficient storage access.
Use a real accessible bucket or create the requested bucket separately through
the user's cloud workflow before retrying. Never silently substitute another
bucket, prefix, or credential. A link's successful creation alone does not prove
that the intended files are readable; use a discovered browse operation to
verify the relevant path when the task requires it.

## Multiple links and result reporting

Prepare the requested name/URI/provider/credential mapping first, then create
links individually through the discovered API. Record returned IDs and statuses
so a partial failure can resume without duplicating successful links. Report
which links succeeded and the concrete error for any failure.

After a mutation, verify it with the discovered lookup or listing operation.
Summarize the workspace, link name, storage URI, status, and requested change.
Include a Platform URL only when it is returned by the API or verified from
the configured Platform context. Do not include authentication or credential
contents in the response.
