#!/bin/bash
echo "../src/SeSiMCMC -i EC-ArcA.dnc -v 0.5 -tr --spaced --html -l+ 22 -r- > arca0501.html 2>arca0501.log"
../src/SeSiMCMC -i EC-ArcA.dnc -v 0.5 -tr --spaced --html -l+ 22 -r- > arca0501.html 2>arca0501.log
echo "../src/SeSiMCMC -i EC-ArcA.dnc -v 0.5 -tr --spaced  --html -l+ 22 -r- -p 0.5 > arca0505.html 2>arca0505.log"
../src/SeSiMCMC -i EC-ArcA.dnc -v 0.5 -tr --spaced  --html -l+ 22 -r- -p 0.5 > arca0505.html 2>arca0505.log
echo "../src/SeSiMCMC -i EC-ArcA.dnc -v 0.1 -tr --spaced --html -l+ 22 -r- > arca0101.html 2>arca0101.log" 
../src/SeSiMCMC -i EC-ArcA.dnc -v 0.1 -tr --spaced --html -l+ 22 -r- > arca0101.html 2>arca0101.log 
echo "../src/SeSiMCMC -i EC-ArcA.dnc -v 0.1 -tr --spaced  --html -l+ 22 -r- -p 0.5 > arca0105.html 2>arca0105.log"
../src/SeSiMCMC -i EC-ArcA.dnc -v 0.1 -tr --spaced  --html -l+ 22 -r- -p 0.5 > arca0105.html 2>arca0105.log


