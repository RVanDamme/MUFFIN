<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# AWS Batch Compute Environment Setup

Complete procedure for setting up an AWS Batch compute environment on Nextflow Platform.

## Prerequisites

- AWS account with appropriate permissions
- AWS CLI configured (`aws sts get-caller-identity` succeeds)
- S3 bucket for Nextflow work directory

## 1. Verify AWS Identity and Permissions

```bash
# Check current identity
aws sts get-caller-identity

# Verify Batch access
aws batch describe-compute-environments --query 'computeEnvironments[].computeEnvironmentName' --output table 2>&1

# Verify S3 access
aws s3api head-bucket --bucket <user-selected-bucket>

# Check current region
aws configure get region
```

### Required IAM Permissions

For Batch Forge (automated provisioning), the IAM user/role needs:

| Service | Permissions | Purpose |
|---------|-------------|---------|
| AWS Batch | Full access (create/describe/update/delete) | Manage compute envs, queues, jobs |
| EC2 | Launch templates, describe VPCs/subnets/SGs | Forge instance provisioning |
| IAM | PassRole, create instance profiles | Delegate roles to Batch instances |
| S3 | Read/write on work dir bucket | Pipeline work directory |
| CloudWatch Logs | Create/describe log groups | Job log streaming |
| ECS | Describe/register task definitions | Container orchestration |

**Optional permissions:**
- EFS / FSx — if using file system storage
- Secrets Manager — if using pipeline secrets
- SSM — for AMI discovery

```bash
# Quick permissions check
aws batch describe-compute-environments 2>&1 | head -3
aws ec2 describe-vpcs --query 'Vpcs[0].VpcId' --output text
aws s3api head-bucket --bucket <work-dir-bucket> 2>&1
```

## 2. Choose Mode: Forge vs Manual

### Batch Forge (Recommended)

Seqera automatically creates Batch compute environments, job queues, IAM roles, and launch templates. Resources are named `TowerForge-*` and can be auto-deleted with "Dispose Resources".

**What Forge creates:**
- Batch compute environment (Spot or On-Demand)
- Job queues (head + compute, separate On-Demand head queue for Spot)
- EC2 launch template with optimized AMI
- IAM instance profile and execution role
- Optional: EFS/FSx file systems, security groups

**Forge only needs:**
- AWS credentials (access key + secret key, or IAM role)
- Region
- S3 work directory
- (Optional) VPC, subnets, instance types

### Manual Mode

User pre-creates all Batch infrastructure and provides queue ARNs.

```bash
# List existing job queues
aws batch describe-job-queues \
  --query 'jobQueues[].{Name:jobQueueName,State:state,Status:status}' --output table

# List existing compute environments
aws batch describe-compute-environments \
  --query 'computeEnvironments[].{Name:computeEnvironmentName,State:state,Status:status,Type:type}' --output table
```

## 3. Configure Work Storage

Use the bucket and work-storage reference selected by the user. Verify access
and region consistency, and create storage only when the user has requested it.
Validate write access with an approved temporary object and remove only that
test object afterward. If lifecycle cleanup is desired, limit its scope to the
user-confirmed intermediate data and retention period; do not prescribe a bucket
name or prefix.

## 4. VPC and Networking (Optional for Forge)

Forge uses default VPC if not specified. For custom networking:

```bash
# List VPCs
aws ec2 describe-vpcs \
  --query 'Vpcs[].{ID:VpcId,CIDR:CidrBlock,Default:IsDefault,Name:Tags[?Key==`Name`].Value|[0]}' \
  --output table

# List subnets in a VPC
aws ec2 describe-subnets \
  --filters "Name=vpc-id,Values=<vpc_id>" \
  --query 'Subnets[].{ID:SubnetId,AZ:AvailabilityZone,CIDR:CidrBlock,Name:Tags[?Key==`Name`].Value|[0]}' \
  --output table

# List security groups
aws ec2 describe-security-groups \
  --filters "Name=vpc-id,Values=<vpc_id>" \
  --query 'SecurityGroups[].{ID:GroupId,Name:GroupName}' \
  --output table
```

### EFS Security Group Requirements (if using EFS)

```bash
# Create security group for EFS
aws ec2 create-security-group \
  --group-name nextflow-efs-sg \
  --description "NFS access for Nextflow" \
  --vpc-id <vpc_id>

# Allow inbound NFS (port 2049) from itself
aws ec2 authorize-security-group-ingress \
  --group-id <sg_id> \
  --protocol tcp --port 2049 \
  --source-group <sg_id>
```

## 5. Create Credentials in Seqera

### Access Keys

```
platform_create_credentials(
    workspace_id=<ws_id>,
    name="aws-<account-alias>",
    provider="aws",
    keys={
        "accessKey": "<aws_access_key_id>",
        "secretKey": "<aws_secret_access_key>"
    }
)
```

### IAM Role (Assume Role)

```
platform_create_credentials(
    workspace_id=<ws_id>,
    name="aws-<account-alias>-role",
    provider="aws",
    keys={
        "assumeRoleArn": "arn:aws:iam::<account_id>:role/<role_name>"
    }
)
```

## 6. Create Compute Environment in Seqera

### Forge Mode

```
platform_create_compute_env(
    workspace_id=<ws_id>,
    name="aws-batch-forge",
    platform="aws-batch",
    credentials_id="<cred_id>",
    config={
        "region": "us-east-1",
        "workDir": "<user-provided work-storage reference>",
        "forge": {
            "type": "SPOT",
            "minCpus": 0,
            "maxCpus": 500,
            "instanceTypes": ["optimal"],
            "disposeOnDeletion": true
        },
        "waveEnabled": true,
        "fusion2Enabled": true
    }
)
```

### Manual Mode

```
platform_create_compute_env(
    workspace_id=<ws_id>,
    name="aws-batch-manual",
    platform="aws-batch",
    credentials_id="<cred_id>",
    config={
        "region": "us-east-1",
        "computeQueue": "arn:aws:batch:us-east-1:<account>:job-queue/<queue>",
        "headQueue": "arn:aws:batch:us-east-1:<account>:job-queue/<head-queue>",
        "workDir": "<user-provided work-storage reference>"
    }
)
```

## Key Configuration Options

| Setting | Default | Notes |
|---------|---------|-------|
| `forge.type` | — | `SPOT` (cheaper) or `EC2` (On-Demand, reliable) |
| `forge.minCpus` | 0 | >0 means always-on instances (costs $) |
| `forge.maxCpus` | — | Upper bound for auto-scaling |
| `forge.instanceTypes` | `["optimal"]` | Or specific types like `["m5.xlarge", "c5.xlarge"]` |
| `waveEnabled` | false | Wave containers for optimized image delivery |
| `fusion2Enabled` | false | Fusion v2 virtual filesystem (requires Wave + NVMe instances) |
| `forge.gpuEnabled` | false | Enable GPU instance types |
| `forge.ebsAutoScale` | true | Auto-expand EBS volumes |

## Troubleshooting

### Credential Issues

```bash
# Verify access key is valid
aws sts get-caller-identity

# Check key permissions
aws batch describe-compute-environments 2>&1

# If using assume role, verify trust policy
aws sts assume-role --role-arn <arn> --role-session-name test 2>&1
```

### Forge Failures

```bash
# Check Batch service-linked role exists
aws iam get-role --role-name AWSServiceRoleForBatch 2>&1
# If not: aws iam create-service-linked-role --aws-service-name batch.amazonaws.com

# Check Batch quota
aws service-quotas get-service-quota \
  --service-code batch \
  --quota-code L-144E0A63 \
  --query 'Quota.Value'
# Default: 50 job queues per account

# Check EC2 instance quotas (Spot)
aws service-quotas get-service-quota \
  --service-code ec2 \
  --quota-code L-34B43A08 \
  --query 'Quota.Value'
```

### S3 Issues

```bash
# Check bucket policy
aws s3api get-bucket-policy --bucket <bucket> 2>&1

# Check bucket encryption (may affect access)
aws s3api get-bucket-encryption --bucket <bucket> 2>&1

# Verify write access with a user-approved temporary object using the selected
# bucket and reference; remove only the object created for that check.
```
