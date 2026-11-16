---
name: nf-plugin-development
description: >
  Create Nextflow plugins using the official nf-plugin-template and publish to the
  Nextflow plugin registry. Covers extension points (operators, functions, factories,
  executors, trace observers), project structure, testing, and publishing. Use when
  the user wants to create a new plugin or add custom functionality to Nextflow.
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


# Nextflow Plugin Development

Create and publish Nextflow plugins using the modern template and registry approach.

## When to Use

Load this skill when the user wants to:
- Create a new Nextflow plugin
- Implement custom executors
- Create trace observers for workflow events
- Publish a plugin to the Nextflow plugin registry
- Test plugins locally before publishing

## Environment Setup

### Verifying Your JDK

Before building, confirm Gradle can see a full JDK (not just a JRE). This is a
common issue in Docker/devcontainer environments where OpenJDK is installed via
`apt` — the `release` file often omits `IMAGE_TYPE=JDK`, causing Gradle's
toolchain resolution to fail with "No matching toolchain found" even though
`java -version` works fine.

Inspect the selected JDK's runtime metadata and confirm Gradle recognizes it as a
full JDK. If it is misclassified, repair or replace that installation using the
appropriate package-management procedure for the user's environment.

When needed, select the supported installed JDK with `JAVA_HOME` before invoking
the project's available Gradle build tool.

## Quick Start: Create a New Plugin

The simplest way to create a plugin:

```bash
# Create a new plugin project (interactive)
nextflow plugin create

# Or specify the name directly
nextflow plugin create nf-my-plugin
```

This scaffolds a complete project based on nf-plugin-template with:
- Gradle build configuration
- Source directories for Groovy/Java code
- Example extension points
- GitHub Actions workflows
- Test infrastructure

## Project Structure

After creation, your plugin has this structure:

Organize the entry workflow, reusable definitions, configuration, helper commands,
tests, and supporting data according to the user's existing project. Locate each
component from the files the user provides rather than assuming a directory layout.

## Extension Points

Plugins can provide several types of extensions. Each requires specific annotations
and interfaces.

### 1. Custom Functions

Functions can be called directly in pipeline scripts:

```groovy
// In MyPluginExtension.groovy
import nextflow.plugin.extension.Function
import nextflow.plugin.extension.PluginExtensionPoint

class MyPluginExtension extends PluginExtensionPoint {

    @Function
    String sayHello(String name) {
        return "Hello, ${name}!"
    }

    @Function
    int addNumbers(int a, int b) {
        return a + b
    }
}
```

Usage in pipeline:

```nextflow
include { sayHello; addNumbers } from 'plugin/nf-my-plugin'

workflow {
    println sayHello('World')      // "Hello, World!"
    println addNumbers(2, 3)       // 5
}
```

### 4. Trace Observers

React to workflow lifecycle events.

**Note:** Use `TraceObserverV2` and `TraceObserverFactoryV2` — the v1 interfaces
are deprecated.

```groovy
import nextflow.trace.TraceObserverV2
import nextflow.trace.TraceObserverFactoryV2
import nextflow.script.WorkflowMetadata

class MyPluginObserver implements TraceObserverV2 {

    @Override
    void onFlowCreate(Session session) {
        // Called from main thread when session starts
        println "Workflow starting: ${session.runName}"
    }

    @Override
    void onFlowComplete(WorkflowMetadata meta) {
        // NOTE: takes WorkflowMetadata parameter — NOT a no-arg method
        // All run-level info (launchDir, runName, success, duration) is on `meta`
        println "Run ${meta.runName} completed. Success: ${meta.success}"
    }

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        // ⚠️ Called from task monitor threads — use thread-safe collections!
        println "Task ${trace.name} completed in ${trace.duration}ms"
    }

    @Override
    void onFilePublish(Path destination, Path source) {
        println "Published: ${destination}"
    }
}

class MyPluginObserverFactory implements TraceObserverFactoryV2 {
    @Override
    Collection<TraceObserverV2> create(Session session) {
        return [new MyPluginObserver()]
    }
}
```

#### `WorkflowMetadata` — What Lives Where

In `TraceObserverV2`, run-level metadata is passed via `WorkflowMetadata meta`
in `onFlowComplete`. Do **not** access these properties on `session` — they
will cause compile errors with `@CompileStatic`:

| Property | Wrong ❌ | Correct ✅ |
|---|---|---|
| Launch directory | `session.launchDir` | `meta.launchDir` |
| Script name | `meta.mainScript` | `meta.scriptName` |
| Run name | `session.runName` | `meta.runName` |
| Revision | `session.revision` | `meta.revision` |
| Commit ID | `session.commitId` | `meta.commitId` |
| Success flag | `session.success` | `meta.success` |
| Duration | — | `meta.duration` |

#### ⚠️ Thread Safety in Observer Callbacks

`onProcessComplete`, `onProcessCached`, and `onProcessSubmit` are called from
**task monitor threads concurrently**, not the main thread. Any mutable state
accumulated in these methods (lists, maps, sets, counters) **must** use
thread-safe types — plain Groovy collections will cause `ConcurrentModificationException`
at runtime, often intermittently and hard to reproduce.

| Use case | Use this |
|---|---|
| Counters | `AtomicInteger` |
| Lists | `new CopyOnWriteArrayList<>()` |
| Maps | `new ConcurrentHashMap<>()` |
| Sets | `new CopyOnWriteArraySet<>()` |

`onFlowCreate` and `onFlowComplete` are called from the main thread and are safe
to use with ordinary variables.

```groovy
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.CopyOnWriteArrayList
import java.util.concurrent.CopyOnWriteArraySet
import java.util.concurrent.atomic.AtomicInteger

class MyPluginObserver implements TraceObserverV2 {

    // Thread-safe state for use in onProcessComplete
    private final AtomicInteger              totalTasks = new AtomicInteger(0)
    private final CopyOnWriteArrayList<String> taskNames = new CopyOnWriteArrayList<>()
    private final ConcurrentHashMap<String, Integer> countByProcess = new ConcurrentHashMap<>()

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        totalTasks.incrementAndGet()
        taskNames.add(trace.name)
        countByProcess.merge(trace.process, 1, Integer::sum)
    }

    @Override
    void onFlowComplete(WorkflowMetadata meta) {
        // Safe to read thread-safe collections here on the main thread
        println "Total tasks: ${totalTasks.get()}"
    }
}
```

### 5. Custom Executors

Implement custom compute backends:

```groovy
import nextflow.executor.Executor
import nextflow.processor.TaskHandler
import org.pf4j.ExtensionPoint

class MyExecutor extends Executor implements ExtensionPoint {

    @Override
    String getName() { return 'my-executor' }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        return new MyTaskHandler(task, this)
    }
}
```

Usage in pipeline:

```nextflow
process MY_PROCESS {
    executor 'my-executor'
    // ...
}
```

## Configuration Scope

Allow users to configure your plugin via `nextflow.config`:

```groovy
import nextflow.plugin.extension.PluginExtensionPoint

class MyPluginConfig {
    String apiUrl = 'https://default.api.com'
    int timeout = 30
    boolean verbose = false
}

class MyPluginExtension extends PluginExtensionPoint {

    @Override
    void init(Session session) {
        def config = session.config.navigate('myplugin') as Map ?: [:]
        this.settings = new MyPluginConfig(config)
    }
}
```

User configuration:

```groovy
// nextflow.config
myplugin {
    apiUrl = 'https://custom.api.com'
    timeout = 60
    verbose = true
}
```

## Registering Extensions


## Building and Testing

### Recommended Development Loop

Follow this order to catch errors at the right layer and avoid wasting time
debugging Nextflow runs when the problem is actually a compile error:

1. Run the project's `compileGroovy` build task to catch syntax and static-type errors.
2. Run its `test` task.
3. Run its local plugin installation task.
4. Run the user-selected validation pipeline with the built plugin.
5. Inspect that run's Nextflow log for silent observer errors.

### Build the Plugin

Use the project's available build interface to run the `assemble` task.

### Test Locally

Install the built plugin using the project's supported local installation task,
then run the user-selected test pipeline with the matching plugin identifier/version.

### Unit Tests

Run the project's `test` build task.

Example test:

```groovy
class MyPluginExtensionTest extends Specification {

    def 'should say hello'() {
        given:
        def ext = new MyPluginExtension()

        expect:
        ext.sayHello('World') == 'Hello, World!'
    }
}
```

### Test Without Publishing

Use `NXF_PLUGINS_TEST_REPOSITORY` to test before publishing:

Set `NXF_PLUGINS_TEST_REPOSITORY` to the repository supplied by the user's build
process, then run the selected validation pipeline with the built plugin version.

## Publishing to the Registry

### 1. Get Registry Access

The user needs to provide an API token to the registry, which they can obtain here: https://registry.nextflow.io/

### 2. Configure Credentials

Configure the `npr.apiKey` Gradle property through the user's supported credential
mechanism. Keep the token out of source files and examples:

```properties
npr.apiKey=your-api-token-here
```

### 3. Prepare for Release

Ensure your project includes:
- A `README.md` file (used as plugin description in registry)
- Correct version in `gradle.properties`
- All tests passing

### 4. Publish

When publishing is authorized, run the project's release or `publishPlugin` task.

The plugin will be available at [registry.nextflow.io](https://registry.nextflow.io).

## Best Practices

1. **Start Simple** — Begin with functions before moving to operators
2. **Follow Naming** — Use `nf-` prefix for plugin names
3. **Test Thoroughly** — Include both unit tests and validation pipeline
4. **Document Well** — Write clear README with usage examples
5. **Version Semantically** — Follow semver for version numbers
6. **Handle Errors** — Provide clear error messages for misconfigurations

## Common Patterns

### Accessing Session Context

```groovy
class MyPluginExtension extends PluginExtensionPoint {

    private Session session

    @Override
    void init(Session session) {
        this.session = session
    }

    @Function
    String getWorkDir() {
        return session.workDir.toString()
    }
}
```

### Working with Files

```groovy
@Function
Path resolveFile(String path) {
    return session.baseDir.resolve(path)
}
```


## Session API Gotchas

### `session.containerConfig` is a Typed Object, Not a Map

`session.containerConfig` returns a typed config object (`DockerConfig`,
`SingularityConfig`, etc.) that does **not** implement `Map`. Accessing it
with map-style property syntax (e.g. `session.containerConfig?.engine`) causes
an `IncompatibleClassChangeError` at runtime.

```groovy
// WRONG — causes IncompatibleClassChangeError at runtime
def engine = session.containerConfig?.engine

// CORRECT — use the getter explicitly, with a safety net for version changes
def engine = 'none'
try { engine = session.containerConfig?.getEngine() ?: 'none' } catch (Throwable ignored) {}
```

### `@CompileStatic` and Groovy Type Strictness

`@CompileStatic` is recommended for all Groovy classes — it catches API errors at
compile time. Watch out for these common Groovy gotchas under static compilation:

```groovy
// Arithmetic on collection results — sum() returns Object, must cast
def avg = (records.collect { it.duration }.sum() as double) / records.size()

// GStrings are not Strings in typed collections — always call .toString()
myStringList << "result-${name}".toString()
myStringMap[key] = "value-${x}".toString()

// Typed map values from computeIfAbsent — use explicit types, not []
myMap.computeIfAbsent(key) { new CopyOnWriteArrayList<SomeType>() }
// NOT: myMap.computeIfAbsent(key) { [] }  — type inference fails
```

## Troubleshooting

### Observer Crashes Are Silent

If your observer throws an exception in `onFlowComplete` or `onProcessComplete`,
Nextflow will **not** fail the pipeline or print the error to stdout. The run
appears to succeed normally — the only sign of failure is that your plugin's
expected output (a file, a notification, etc.) is missing.

**Always check `.nextflow.log` after a test run:**

```bash
grep -i "exception\|error" .nextflow.log | grep -v "DEBUG"
```

To make errors visible during development, wrap your `onFlowComplete` body:

```groovy
@Override
void onFlowComplete(WorkflowMetadata meta) {
    try {
        writeOutput(meta)
    } catch (Throwable e) {
        log.error("nf-my-plugin: failed to write output — ${e.message}", e)
        // Optionally rethrow to surface the error in the pipeline exit code:
        // throw new RuntimeException("nf-my-plugin failed", e)
    }
}
```

### Plugin Not Loading

- Verify `build.gradle` specifies the plugin entrypoint and extension classes
- Ensure version matches in all config files

### Functions Not Found

- Confirm `@Function` annotation is present
- Check the function is public
- Verify plugin is included: `-plugins nf-my-plugin`

### Tests Failing

- Run with `--stacktrace` for detailed errors
- Check Nextflow version compatibility

## Related Skills

- `nf-plugin-legacy-migration` — Convert legacy plugins to new registry format
- `nextflow-config` — Configure plugins in nextflow.config
- `nf-pipeline-design` — Standard pipeline organization
