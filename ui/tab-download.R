tabPanel("Download",
         fluidRow(
           box(
             title = "Training Set",
             width = 6, background = "light-blue",
             div(downloadButton(outputId = "dlf1",
                                label = "Download",
                                href = "https://drive.google.com/uc?export=download&id=0B9o4VIJJSodfMmd0Um9LWXM3VFE",
                                icon = icon("download")),
                 align = "center")
           ),
           box(
             title = "Evaluation Set",
             width = 6, background = "yellow",
             div(downloadButton(outputId = "dlf2",
                                label = "Download",
                                href = "https://drive.google.com/uc?export=download&id=0B9o4VIJJSodfNlZMNHpPWW1oSUE",
                                icon = icon("download")),
                 align = "center")
           )
         ),
         fluidRow(
           box(
             title = "Real dataset",
             width = 6, background = "red",
             div(downloadButton(outputId = "dlf3",
                                label = "Download",
                                href = "https://drive.google.com/uc?export=download&id=0B9o4VIJJSodfOTNOelRKNm9MQ2M",
                                icon = icon("download")),
                 align = "center")
           )
         )
)