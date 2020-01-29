#!/usr/bin/env bash

clean_tcoffee(){
    (find /tmp -name "t_coffee.tmp" -print0 \
    | xargs dirname \
    | xargs rm -rf \
    ) || true
}


clean_tcoffee