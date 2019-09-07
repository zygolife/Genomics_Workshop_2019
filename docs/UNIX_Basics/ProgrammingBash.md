# Shell programming in BASH

The shell is powerful environment and can be used to encode lots of logic.  At some point it can be useful to switch to even more expressive programming and scripting language like 


# Logical operators

**if / elif / else**

```
if [ ! -f file.txt ]; then
    echo "Gene X" > file.txt
fi
```

# iterator

```
for file in $(ls *.fa)
do
  base=$(basename $file .fa) # strip the .fa from the end of the filename
  blastp -query $file -db swissprot -out $base.BLASTP.tab -outfmt 6
done
```


