# Stage logs

One file per stage (`stage-<N>.md`), append-only. One entry per task:

```
## <task id> — <title> — <YYYY-MM-DD HH:MM>
Changed: <files or summary>
Checks:
  - <command> → <result, with the relevant output lines>
Not verified: <what and why, or "nothing">
Decisions: <any, also copied to PROGRESS.md>
Previous entry commit: <hash>
```

Saved stage reports go in `../reports/stage-<N>.md`.
