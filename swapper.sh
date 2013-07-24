#!/bin/bash
# Manually create swap file

SWAP_SIZE_MEGABYTES=2000
    if [ $SWAP_SIZE_MEGABYTES -eq 0 ];then
      echo No swap size given, skipping.
   else
     if [ -e /mnt/swapfile ];then
        echo /mnt/swapfile already exists.  Skipping.
     else
        echo Creating /mnt/swapfile of $SWAP_SIZE_MEGABYTES Megabytes
        dd if=/dev/zero of=/mnt/swapfile bs=1024 count=$(($SWAP_SIZE_MEGABYTES*1024))
        mkswap /mnt/swapfile
        swapon /mnt/swapfile
        echo Swap Status:
        swapon -s
    fi
  fi

