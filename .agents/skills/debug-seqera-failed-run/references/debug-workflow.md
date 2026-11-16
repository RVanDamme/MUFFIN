<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Nextflow Platform Run Debug Workflow

Use this workflow when diagnosing a failed Nextflow Platform run.

## 1. Resolve the workflow

```text
platform_list_workflows(workspace_id=<id>, max=1)
```

If the user specifies a run name or workflow ID, use that instead. If no
workspace is known, call `platform_list_user_workspaces` first.

## 2. Get workflow details and progress

```text
platform_get_workflow(workflow_id=<id>, workspace_id=<id>)
platform_get_workflow_progress(workflow_id=<id>, workspace_id=<id>)
```

Check `status`, `exitStatus`, `errorMessage`, `start`, `complete`, `duration`,
and task counts.

## 3. Find failed or incomplete tasks

```text
platform_list_workflow_tasks(workflow_id=<id>, workspace_id=<id>, status="FAILED")
```

If no explicit status filter is available, list all tasks and filter by
`status != "SUCCEEDED"`.

## 4. Inspect each failed task

```text
platform_get_workflow_task(workflow_id=<id>, task_id=<id>, workspace_id=<id>)
platform_get_workflow_task_log(workflow_id=<id>, task_id=<id>, workspace_id=<id>)
```

Check `exitStatus`, `stderr`, `stdout`, `process`, `tag`, `cpus`, `memory`,
`realtime`, `pcpu`, `pmem`, `rchar`, and `wchar`.

## 5. Read workflow-level logs when task logs are insufficient

```text
platform_get_workflow_log(workflow_id=<id>, workspace_id=<id>)
```

Workflow-level logs are especially important for failures before task execution:
script compilation, config validation, malformed sample sheets, credential
errors, plugin/download errors, channel errors, and staging failures.
