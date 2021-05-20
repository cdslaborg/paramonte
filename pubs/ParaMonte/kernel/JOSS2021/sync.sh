#!/bin/bash
manFiles="paper.bib:paper.md:PredictedSpeedupActualSpeedup.png:procContribution512.png"

for file in ${manFiles//:/ }; do
    filePath="../../../../../manParaMonteJOSS2020/current_JOSS/${file}"
    echo >&2
    echo >&2 "-- copying manuscript file,"
    echo >&2 "-- from: ${filePath}"
    echo >&2 "-- to: ./"
    echo >&2
    cp "${filePath}" "./"
done
