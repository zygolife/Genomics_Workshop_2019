# Running programs


# Redirect output

To save the output of a program use `>` which will create a new file with this name. If this file already exists **THIS WILL OVERWRITE THE EXISTING FILE**.
```
echo "This is my song" > song.txt
```
If you want to **APPEND** to the end of a file use ``>>``.
```
echo "After killing" > haiku.txt
echo "a spider, how lonely I feel" >> haiku.txt
echo "in the cold of night!" >> haiku.txt
```


# Pipe commands to connect applications together

The ouput of one program can be passed as the input to another program.

```command1 | command2```
The output from command1 goes to command 2

# Shell commands and useful keys

  * `^` means "Control key"
  * cancel a running application: `^C`
  * end a session: `^D` (End of File message)
  * Push the `Tab` key to try to autocomplete (applications, filenames, directories)
  * while typing on command line - jump to end of line: `^E`
  * get back to the beginning of line: `^A`
  * Up and down keys cycle through history of commands
  * Type `!!` to execute the last command
  * Type `history` to see list of previous commands
  * Type !NUMBER to excute cmd from that list
