graph [
  directed 1
  node [
    id 0
    label "VHL"
    class "driver"
  ]
  node [
    id 1
    label "HIF1A"
    class "other"
  ]
  node [
    id 2
    label "EPAS1"
    class "other"
  ]
  node [
    id 3
    label "SLC2A1"
    class "effector"
  ]
  node [
    id 4
    label "VEGFA"
    class "effector"
  ]
  node [
    id 5
    label "PDGFB"
    class "effector"
  ]
  node [
    id 6
    label "TGFA"
    class "effector"
  ]
  node [
    id 7
    label "CCND1"
    class "effector"
  ]
  edge [
    source 0
    target 1
    effect "inhibition"
  ]
  edge [
    source 0
    target 2
    effect "inhibition"
  ]
  edge [
    source 1
    target 3
    effect "stimulation"
  ]
  edge [
    source 1
    target 4
    effect "stimulation"
  ]
  edge [
    source 1
    target 5
    effect "stimulation"
  ]
  edge [
    source 1
    target 6
    effect "stimulation"
  ]
  edge [
    source 2
    target 7
    effect "stimulation"
  ]
]
