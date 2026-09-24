#import "@preview/quill:0.8.0": *

#set page(width: auto, height: auto, margin: 10pt)
#set text(font: "New Computer Modern", size: 11pt)

#align(center)[
  #quantum-circuit(
    lstick($p_1':|0〉$), $H$, ctrl(1), ctrl(2),[\ ],
    lstick($p_1'':|0〉$),[\ ],
    lstick($p_2':|0〉$),[\ ],
    lstick($p_2'':|0〉$),[\ ],
    lstick($p_3':|0〉$),[\ ],
    lstick($p_4':|0〉$),[\ ],
    lstick($p_4'':|0〉$),[\ ],  
    lstick($p_5':|0〉$),[\ ],
    lstick($p_5'':|0〉$),[\ ],
    lstick($p_6':|0〉$),[\ ],
    lstick($p_6'':|0〉$),[\ ],
    lstick($p_7':|0〉$),[\ ],
    lstick($p_8':|0〉$),[\ ],
    lstick($p_8'':|0〉$),[\ ],
    lstick($p_9':|0〉$),[\ ],
    lstick($p_9'':|0〉$),[\ ],

    lstick($a_(h_1):|0〉$),[\ ],
    lstick($a_(h_2):|0〉$),[\ ],
    lstick($a_(v_1):|0〉$),[\ ],
    lstick($a_(v_2):|0〉$),[\ ],
    lstick($a_(v_3):|0〉$),[\ ],
  )
]