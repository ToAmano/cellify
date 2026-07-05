---
name: safe-push
description: コードを整形し、テストを通してから安全にリモートへブランチをプッシュする
---

# 実行手順
1. `pre-commit run --all-files` を実行してコードを整形する。
2. `pytest` でテストを実行する。
3. すべてパスしたら、現在の変更をコミットして `git push origin <ブランチ名>` を実行する。
