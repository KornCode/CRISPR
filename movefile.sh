#!/bin/sh

userName=$(whoami)
masterIP=$(echo $SSH_CONNECTION | awk '{print $1}')

scp reverse.txt $userName@$masterIP:/project

if [ $? -eq 0 ]; then
	echo File from worker has been sent to master successfully.
else
	echo File from worker cannot send to master.
fi


