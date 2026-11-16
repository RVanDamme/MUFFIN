---
name: nf-plugin-legacy-migration
description: >
  Migrate legacy Nextflow plugins to the new plugin registry format. Covers converting
  plugins that use the old GitHub-based nf-plugins-index to the modern registry.nextflow.io
  approach. Use when the user has an existing plugin using the legacy format or wants
  to understand the differences between old and new plugin systems.
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


# Nextflow Plugin Legacy Migration

Migrate a user's legacy Nextflow plugin to the modern registry format without assuming a local project layout.

## When to Use

Use for plugins based on the older nf-hello template or GitHub plugin index, migration to nf-plugin-template, registry publishing, and compatibility across Nextflow versions.

## Assess the Existing Plugin

Inspect the source and build files supplied by the user. Identify the plugin class, extension classes, observer interfaces, package namespace, Gradle configuration, extension registration, release metadata, tests, and current publishing mechanism.

Legacy indicators include nf-hello ancestry, references to the GitHub plugin index, older Gradle plugin syntax, and a release process that updates index metadata manually.

The bundled baseline distinguishes the legacy index used before Nextflow 25.10
from the modern registry approach in 25.10 and later. Verify support for the
user's intended release rather than inferring it from the template alone.

## Migration Workflow

1. Confirm the user's supported Nextflow versions and intended registry package name.
2. Compare the existing project with the current official nf-plugin-template. A fresh template can be useful, but agree how it fits the user's project before moving or replacing files.
3. Preserve the plugin's behavior while migrating source classes and package names.
4. Update deprecated observer interfaces to the interfaces supported by the target Nextflow version. Check both the observer and its factory.
5. Update extension registration to contain the actual fully qualified extension class names.
6. Adapt the user's Gradle configuration, release metadata, and existing build commands to the current template.
7. Preserve and run the existing unit tests; add validation for behavior affected by the migration.
8. Install the migrated plugin in the user's chosen Nextflow environment and test a small workflow.
9. Configure registry authentication through the supported local credential mechanism. Never embed a real API token in source, examples, or responses.
10. Publish only when the user authorizes publishing, then verify the registry entry.

## Source Changes

Use the user's organization namespace consistently in source, tests, and registration metadata. Nextflow extension annotations remain conceptually distinct:

~~~groovy
@Function
String myFunction(String argument) { /* implementation */ }

@Operator
DataflowWriteChannel myOperator(DataflowReadChannel source) { /* implementation */ }

@Factory
DataflowWriteChannel myFactory(Map parameters) { /* implementation */ }
~~~

When migrating observer APIs, compare each implemented method with the current interface instead of changing only the class declaration. Initialization callbacks and file-publication callbacks may differ between versions.

For the V2 observer baseline, migrate TraceObserver and its factory to
TraceObserverV2 and TraceObserverFactoryV2 respectively. Check the callbacks:

| Legacy callback | V2 baseline |
|---|---|
| onFlowStart(Session) | onFlowCreate(Session) |
| onFlowComplete() | onFlowComplete() |
| onProcessComplete(TaskHandler, TraceRecord) | Same |
| onProcessCached(TaskHandler, TraceRecord) | Same |
| onFilePublish(Path) | onFilePublish(Path, Path) |

Configuration access should continue to respect user settings and missing values. For example:

~~~groovy
def config = session.config.navigate('myplugin') as Map ?: [:]
~~~

## Build and Release Configuration

Use the current io.nextflow.nf-plugin-gradle version supported by the user's target Nextflow versions. Keep package name, plugin ID, version, provider, minimum Nextflow version, plugin class, and extension registrations consistent.

When the user's Gradle build needs public plugin repositories, the template can
use the following configuration without prescribing a local project layout:

~~~groovy
pluginManagement {
    repositories {
        gradlePluginPortal()
        maven { url 'https://plugins.gradle.org/m2/' }
    }
}
~~~

The registry description should explain installation, usage, configuration, and licensing. Keep the user's existing build interface where possible; do not assume a particular Makefile or working directory exists.

## Compatibility and Troubleshooting

For older supported Nextflow releases, verify actual API availability. Do not infer compatibility from a version comparison alone.

If a migrated plugin is not found, check its published package name and version, registry configuration, installed metadata, and the user's existing plugin cache. Propose any cache cleanup using the location identified in that environment, with scope limited to the affected plugin.

If observers are not called, inspect interface compatibility and factory registration. If extensions are missing, verify fully qualified registration names. For Gradle failures, compare the wrapper and build plugin requirements before changing versions.

## Completion Checklist

- Source and registration metadata agree.
- Supported Nextflow versions have been checked.
- Existing behavior is preserved.
- Unit tests and local plugin installation pass.
- A representative workflow runs.
- Release metadata and documentation describe the migrated plugin.
- Any publishing action is authorized.

## Public References

- https://github.com/nextflow-io/nf-plugin-template
- https://registry.nextflow.io

Use nf-plugin-development for new plugins and nextflow-config for pipeline plugin configuration.
