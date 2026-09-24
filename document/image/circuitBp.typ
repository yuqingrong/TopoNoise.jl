#import "@preview/quill:0.8.0": *

#set page(width: auto, height: auto, margin: 10pt)
#set text(font: "New Computer Modern", size: 11pt)

#align(center)[
  #quantum-circuit(
    lstick($p_1:|0〉$),1,targ(),  [\ ],
    lstick($p_2:|0〉$),3,targ(),  [\ ],
    lstick($p_3:|0〉$),           [\ ],
    lstick($a_1:|0〉$),$H$, ctrl(-3), ctrl(1), 1,swap(3),5,swap(6),[\ ],
    lstick($a_2:|0〉$),2,targ(),ctrl(-3),1,swap(3), $H$,ctrl(3),ctrl(1),2,swap(6),[\ ],
    lstick($a_3: |0〉$),8,targ(),ctrl(3),2,swap(6),[\ ],
    lstick($p_4: |0〉$),4,swap(),[\ ],
    lstick($p_5: |0〉$),5,swap(),1,targ(),[\ ],
    lstick($p_6: |0〉$),9,targ(),[\ ],
    lstick($p_7: |0〉$),10,swap(),[\ ],
    lstick($p_8: |0〉$),11,swap(),[\ ],
    lstick($p_9: |0〉$),12,swap(),
  )
]
