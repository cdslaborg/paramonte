#!/bin/bash
manFiles="adaptationMeasure.png:covmatEvolution.png:paper.bib:paper.md:PredictedSpeedupActualSpeedup.png:procContribution512.png"

for file in ${manFiles//:/ }; do
    filePath="../../../../../manParaMonteJOSS2020/${file}"
    echo >&2
    echo >&2 "-- copying manuscript file,"
    echo >&2 "-- from: ${filePath}"
    echo >&2 "-- to: ${filePath}"
    echo >&2
    cp "${filePath}" "./"
done
