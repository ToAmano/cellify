# Custom Agent Rules for cellify

This document defines project-scoped rules and style guidelines for AI agents working in this repository.

## 1. Documentation Synchronization

Whenever you modify any features, CLI arguments, core logic, or adapters in `cellify`:
- **README.md**: You MUST update [README.md](file:///Users/amano/works/research/supercell/README.md) to reflect the changes (e.g., new arguments, examples, or behavior changes).
- **SKILL.md**: You MUST update the Antigravity Skill guide at [.agents/skills/cellify/SKILL.md](file:///Users/amano/works/research/supercell/.agents/skills/cellify/SKILL.md) to ensure other agents can correctly utilize the updated features.

## 2. Safe Push Rules

Before committing or pushing any code to the remote repository, you MUST:
1. Run `pre-commit run --all-files` to format and lint the code.
2. Run `pytest` to ensure all unit tests pass successfully.
3. Never commit or push if formatting, linting, or tests are failing.
