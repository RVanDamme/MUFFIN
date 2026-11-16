<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Seqera Compute Environment Setup

Complete procedure for setting up a Seqera Compute environment on Nextflow Platform.

Seqera Compute is a fully managed compute environment available on Seqera Cloud. Seqera automatically provisions and manages all underlying resources including AWS accounts, credentials, roles, compute environments, and S3 storage buckets. Users do not need to provide cloud credentials or configure infrastructure.

## Prerequisites

- A Seqera Cloud account (Seqera Compute is **not available** on self-hosted Platform instances)
- A workspace on Seqera Cloud

**Not required:**
- No AWS account needed
- No cloud credentials needed
- No S3 bucket setup needed
- No VPC/networking configuration needed

## Key Differences from Other Platforms

| Aspect | AWS Batch / Kubernetes | Seqera Compute |
|--------|----------------------|----------------|
| Credentials | User provides cloud credentials | None needed (managed by Seqera) |
| Infrastructure | User manages cloud resources | Fully managed by Seqera |
| Work directory | User provides S3/cloud path | Auto-provisioned S3 bucket |
| Configuration | Many required parameters | Only name and region required |
| Availability | Cloud + self-hosted | Seqera Cloud only |

## Configuration Options

| Field | Required | Description |
|-------|----------|-------------|
| `name` | Yes | Descriptive name (e.g., "Seqera Compute eu-west-1") |
| `region` | Yes | AWS region for execution |
| `workDir` | No | User-selected storage suffix appended to managed storage; use the value type required by the discovered API rather than a full storage reference. |
| `preRunScript` | No | A literal shell command (not a file path) executed before workflow tasks. |
| `postRunScript` | No | A literal shell command executed after workflow tasks. |
| `nextflowConfig` | No | Additional Nextflow configuration |
| `environment` | No | Environment variables. Array of `{name, value, head, compute}`. `head` exposes the var to the workflow controller; `compute` exposes it to task containers. Most users want both `true`. |

Other optional fields (e.g., instance sizing for free-tier vs paid orgs) exist on the Platform side. Confirm the exact field names against the live API rather than guessing — use `search_seqera_api` with a query like "create compute env" to see the current suggested parameters.

### Supported Regions

Use `search_seqera_api` with a query like "list platform regions" to discover the right tool and call it dynamically — the api_names change as the MCP server evolves. Common regions include:

- `us-east-1` (N. Virginia)
- `eu-west-1` (Ireland)
- `eu-west-2` (London)

## Create Compute Environment

### Via MCP Tool

```
platform_create_compute_env(
    workspace_id=<ws_id>,
    name="seqera-compute-eu",
    platform="seqeracompute-platform",
    config={
        "region": "eu-west-1"
    }
)
```

The `platform` value is `"seqeracompute-platform"` — `"seqera-compute"` and `"seqera_compute"` both look correct and both get rejected.

Omit `credentials_id` — Seqera Compute provisions and uses its own credentials internally, so passing one is unnecessary and produces a confusing auth error rather than overriding anything useful.

### With Optional Settings

```
platform_create_compute_env(
    workspace_id=<ws_id>,
    name="seqera-compute-prod",
    platform="seqeracompute-platform",
    config={
        "region": "us-east-1",
        "workDir": "<user-selected storage suffix>",
        "preRunScript": "export NXF_OPTS='-Xms1g -Xmx4g'",
        "environment": [
            {"name": "MY_VAR", "value": "my-value", "head": true, "compute": true}
        ]
    }
)
```

## Common Pitfalls

| Symptom | Cause | Fix |
|---------|-------|-----|
| "Platform not available" | Running against self-hosted Platform | Seqera Compute only works on Seqera Cloud (cloud.seqera.io) |
| Auth/credential error | Passed `credentials_id` | Drop it — Seqera Compute provisions and uses its own credentials internally |
| Platform-ID rejected | Used `"seqera-compute"` or `"seqera_compute"` | The exact ID is `"seqeracompute-platform"` |
| `workDir` validation fails | Value does not match the API contract | Confirm the suffix expected by the discovered managed-compute schema |
| Invalid region | Region not supported on Seqera Compute | Discover supported regions via `search_seqera_api` rather than guessing |
| Free-tier rejection | Instance size above the free-tier cap | Free-tier orgs are limited to the smallest instance size; check current limits against the live API |
