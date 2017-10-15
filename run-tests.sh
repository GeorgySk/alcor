#!/usr/bin/env bash

docker-compose -f docker-compose.db.yml -f docker-compose.yml up --exit-code-from alcor
STATUS=$?

docker-compose down

if [ "$STATUS" -eq "0" ]; then
	echo "Tests passed";
else
	echo "Tests failed to pass"
fi

exit ${STATUS}
