## Interacting with Text Files

UNIX and some of the commong scripting languages (Python, Perl, Ruby) are adept at interacting with and processing text files.

Reminder that  **more** and **less* are ways to read these files one page at a time

**grep**


### Text Editors

**nano** also try **pico** - these are super simple text editors with build in menus. Useful on the command line.

**vi** - One of the most powerful editors. It has a lot of built-in commands to remember which require some mastery. Most common question is how to exit -- use 'ESC' and then type ":q" - the `:` is part of telling vi you are sending it a command. `q` means quite and `wq` means write and quit.

There are many useful built in commands
   - use `/` to search for a term.
   - use `:s/x/y/` to replace all instances of `x` on a line with `y`
   - use `:s/x/y/g` to globally replace all instances of x with y on a line

**emacs** - A powerful graphical or command-line editor. There are graphical menus.

   - Try `emacs -nw` to run emacs without graphical window. Useful on a command-line only login interface.
