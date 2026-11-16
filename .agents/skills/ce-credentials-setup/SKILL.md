---
name: ce-credentials-setup
description: >
  Set up Nextflow Platform compute environments (Seqera Compute, AWS Batch, Kubernetes)
  and the cloud or cluster credentials they need. Walks through gathering prerequisites
  — IAM keys, bearer tokens, certificates, regions, work directories — and then
  calling platform_create_credentials and platform_create_compute_env. Use this skill
  whenever a user wants to wire Nextflow Platform to any compute backend, troubleshoots
  a failing CE, asks about Seqera Compute as a managed option, or talks about
  "connecting my AWS account", "adding my cluster", or "using my own compute" — even
  if they don't explicitly say "compute environment" or "credentials".
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


# Compute Environment & Credentials Setup

Step-by-step procedures for gathering prerequisites before creating CEs and credentials on Nextflow Platform.

## When to Use

Load this skill when the user wants to:
- Set up a new compute environment on Nextflow Platform
- Create cloud or Kubernetes credentials
- Understand what values are needed for a specific platform
- Troubleshoot CE creation failures (wrong token, missing permissions, etc.)

Read the [MCP connection](../launch-workflow/references/seqera-mcp/README.md) connection reference for discovery and authentication. This skill owns provider prerequisites and compute readiness.

## General Workflow

1. **Identify platform** — ask which cloud/cluster (seqera-compute, aws-batch, k8s, etc.)
2. **Gather credentials** — run commands to get keys, tokens, certificates (skip for Seqera Compute)
3. **Verify locally** — confirm credentials work before sending to Platform (skip for Seqera Compute)
4. **Create credentials** — `platform_create_credentials` through discovered MCP operations (skip for Seqera Compute)
5. **Create compute environment** — `platform_create_compute_env` through discovered MCP operations

> **Seqera Compute shortcut:** If the user is on Seqera Cloud and wants the simplest option, recommend Seqera Compute. It only needs a name and region — no credentials, no infrastructure setup. Skip directly to step 5. This is the path with the fewest moving parts and the fewest ways to fail.

## Credential Types Quick Reference

| Provider | Required Values | Verify Command |
|----------|----------------|----------------|
| Seqera Compute | None (fully managed) | N/A |
| AWS | Access Key, Secret Key (or IAM Role ARN) | `aws sts get-caller-identity` |
| Kubernetes | API server URL, Bearer Token, SSL Certificate | `kubectl cluster-info` |
| Container Registry | Registry URL, Username, Password | `docker login <registry>` |
| SSH | Private key, passphrase (optional) | Verify the user-selected host with the approved SSH identity |

## Seqera Compute

Fully managed compute — Seqera provisions all infrastructure automatically. No credentials needed.

**Requirements:** Seqera Cloud account only (not available on self-hosted Platform).

Invoke via `call_seqera_api`:

```json
{
  "service": "platform",
  "api_name": "platform_create_compute_env",
  "parameters": {
    "workspaceId": 12345,
    "name": "seqera-compute-eu",
    "platform": "seqeracompute-platform",
    "config": {"region": "eu-west-1"}
  }
}
```

For full configuration options, instance sizes, and troubleshooting: load `seqera-compute.md` via the host's file-reading tool.

## Kubernetes

Requires: API server URL, service account token, CA certificate, namespace with `ReadWriteMany` PVC.

### Quick Setup

```bash
# 1. Verify cluster access
kubectl cluster-info

# 2. Apply Seqera manifest (creates namespace, service account, roles).
#    Find the current URL at https://github.com/seqeralabs/nf-tower-k8s/releases —
#    do not blindly pin `master`, the file has been renamed historically.
kubectl apply -f <manifest_url_from_latest_release>

# 3. Get bearer token (secret name comes from the manifest used in step 2)
kubectl describe secrets/<token-secret-name> -n <tower-namespace>

# 4. Get API server URL
kubectl cluster-info | grep "Kubernetes control plane"

# 5. Get SSL certificate
kubectl config view --raw -o jsonpath='{.clusters[0].cluster.certificate-authority-data}' | base64 -d
```

For full procedure including RBAC, storage setup, and troubleshooting: load `kubernetes.md` via the host's file-reading tool.

## AWS Batch

Two modes: **Batch Forge** (Seqera creates resources automatically) or **Manual** (user provides existing queues).

### Quick Setup (Forge)

```bash
# 1. Verify AWS access
aws sts get-caller-identity

# 2. Check required permissions exist
aws batch describe-compute-environments 2>&1 | head -5

# 3. Identify target region
aws configure get region

# 4. List available VPCs/subnets (optional — Forge uses defaults)
aws ec2 describe-vpcs --query 'Vpcs[].{ID:VpcId,CIDR:CidrBlock,Default:IsDefault}' --output table
aws ec2 describe-subnets --query 'Subnets[].{ID:SubnetId,VPC:VpcId,AZ:AvailabilityZone}' --output table

# 5. Verify access to the user's selected work-storage bucket
aws s3api head-bucket --bucket <user-selected-bucket>
```

For full procedure including IAM policies, Forge vs Manual (use Forge unless your team requires explicit IAM control over the Batch queues), EFS/FSx, and Fusion: load `aws-batch.md` via the host's file-reading tool.

## Calling MCP Tools

After gathering values, create credentials then CE. Both go through `call_seqera_api` (the authenticated MCP connection handles credentials — never pass auth headers yourself).

Step 1: Create credentials.

```json
{
  "service": "platform",
  "api_name": "platform_create_credentials",
  "parameters": {
    "workspaceId": 12345,
    "name": "my-k8s-creds",
    "provider": "k8s",
    "keys": {
      "token": "<bearer_token>",
      "certificate": "<ca_cert_pem>"
    }
  }
}
```

Step 2: Create compute environment, passing the credentials ID from step 1.

```json
{
  "service": "platform",
  "api_name": "platform_create_compute_env",
  "parameters": {
    "workspaceId": 12345,
    "name": "k8s-prod",
    "platform": "k8s",
    "credentialsId": "<cred_id_from_step_1>",
    "config": {
      "server": "<api_server_url>",
      "namespace": "tower-nf",
      "headServiceAccount": "tower-launcher-sa",
      "storageClaimName": "tower-scratch",
      "storageMountPath": "<user-provided mount location>",
      "workDir": "<user-provided work location>"
    }
  }
}
```

Use `search_seqera_api` to confirm exact field names against the live API before calling — Platform's REST schema evolves and the MCP layer's normalization can drift from what's documented here.

## Troubleshooting

| Error | Cause | Fix |
|-------|-------|-----|
| Invalid credentials | Wrong token or expired | Re-run token extraction command, check expiry |
| Permission denied | Missing IAM policy or RBAC role | Check required permissions in reference doc |
| Region mismatch | Credential region ≠ CE region | Ensure consistent regions in both |
| SSL certificate error | Self-signed, expired, or wrong cert | Re-extract CA cert from kubeconfig or cluster |
| Connection refused | Wrong API server URL or firewall | Verify URL with `kubectl cluster-info`, check network |
| Storage not found | PVC doesn't exist or wrong name | Verify PVC with `kubectl get pvc -n <namespace>` |
| Quota exceeded | AWS Batch queue limit (50/account) | Delete unused CEs or request quota increase |
