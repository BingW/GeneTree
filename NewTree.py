from ete2 import Tree, TreeStyle

t = Tree("(((((((((((Scer|YLL043W,Spar|spar53-g24.1),Smik|smik1363-g1.1),Sbay|sbayc605-g4.1),Cgla|CAGL0C03267g),Scas|Scas400.1),(((((Scer|YFL054C,Spar|spar345-g20.1),Smik|smik1154-g4.1),Sbay|sbayc622-g21.1),Cgla|CAGL0E03894g),Scas|Scas622.5)),((Kwal|Kwal55.20572,((Klac|KLLA0E00550g,Sklu|SAKL0H25916g),Agos|ACL068W)),Kwal|Kwal33.15269)),Clus|CLUG05700),(Ylip|YALI0F00462g,Ylip|YALI0E05665g)),(Anid|AN3915,Anid|AN0830)),Spom|SPAC977.17);")
circular_style = TreeStyle()
circular_style.mode = "c"
circular_style.scale = 20
t.render("/Users/bingwang/VimWork/GeneTree/c_test.png",w=183,units="mm",tree_style=circular_style)

