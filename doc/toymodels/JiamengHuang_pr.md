A new parameter "vol_abs" is added. If you want to use absolute volume to get EOS, you can add

         "vol_abs":      true,

in the "eos" part of property.json
if it's not mentioned, "False" is set defaultly
when you are using absolute volume, there will be a notation in the last line of output during "make" process, which is like

(absolute volume)
