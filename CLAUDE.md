# Stats::LikeR — instructions for Claude

## Never edit `Changes`

Do not create, edit, append to, or revert `Changes` under any circumstances.
The release notes there are written by hand, by the maintainer, and are not
Claude's to touch — not even to add an entry for work Claude just did, and not
even when asked to "update the changelog" as part of a larger task.

This holds regardless of the tool: no `Edit`, no `Write`, no `sed -i`/`perl -pi`,
no `git checkout`/`git revert` that touches it, no patch that includes it.

When a change would normally warrant a release note, say so in the reply and
leave the wording to the maintainer. Do not draft it into the file.
