# Interacting with Text Files

UNIX and some of the commong scripting languages (Python, Perl, Ruby) are adept at interacting with and processing text files.

Reminder that  **more** and **less* are ways to read these files one page at a time

**grep** allows for searching for a text string in one or multiple files. There are powerful expressions that can be used.

`-A[N]` - display N lines after a match
`-B[N]` - display N lines before a match
`-P`    - support powerful perl / python regular expression syntax
`-c`    - return the count of the number of lines that match
`-v`    - only report lines that don't match

```bash
grep -c "\texon\t" data/genomic/Ncrassa.gff
grep -c "\tmRNA\t" data/genomic/Ncrassa.gff
10815
```

# Text Editors

**nano** also try **pico** - these are super simple text editors with build in menus. Useful on the command line.

**vi** - One of the most powerful editors. It has a lot of built-in commands to remember which require some mastery. Most common question is how to exit -- use 'ESC' and then type ":q" - the `:` is part of telling vi you are sending it a command. `q` means quite and `wq` means write and quit.

There are many useful built in commands
   - use `/` to search for a term.
   - use `:s/x/y/` to replace all instances of `x` on a line with `y`
   - use `:s/x/y/g` to globally replace all instances of x with y on a line

**emacs** - A powerful graphical or command-line editor. There are graphical menus.

   - Try `emacs -nw` to run emacs without graphical window. Useful on a command-line only login interface.
   - Graphical windows make it easy to switch back and forth between command line.
   - When you start up graphical emacs put it in the background by typing the `&` at end.
  So
  ```Bash
  emacs myfile.sh &
  ```

# Formatted printing

On bash can can use the **echo** and **printf** commands to print messages. These can be sent to the screen or to the

**printf** can be used
```
printf "%-10s hello\n" "ab" # left formatted
printf "%10s hello\n" "ab" # right formatted
printf "%.2f formatted float\n" "12.333334"
printf "%.2f\tformatted float\n" "12.333334"

$ printf "%.2e formatted scientific\n" "120000.1"
1.20e+05 formatted scientific
```
