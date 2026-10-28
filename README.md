# MammoStationarity

**Aplicativo em Python para análise da distribuição espacial do ruído em imagens médicas.**

O MammoStationarity é uma ferramenta computacional desenvolvida para investigar a variabilidade do ruído em imagens médicas por meio de métricas estatísticas e visualização de dados.

O aplicativo realiza uma varredura sistemática da imagem utilizando regiões de interesse (ROIs), calcula estatísticas locais e apresenta os resultados em mapas de calor e gráficos. A ferramenta foi desenvolvida no contexto de pesquisa em Física Médica, com foco inicial em mamografia digital.

## Funcionalidades

- **Leitura de imagens DICOM:** carregamento dos dados de pixel e utilização das informações de espaçamento e transformação de valores disponíveis nos metadados.
- **Varredura bidimensional:** análise de regiões de interesse distribuídas ao longo da imagem.
- **Análise estatística local:** cálculo do desvio padrão e da média dos valores dos pixels em cada ROI.
- **Mapas de calor:** visualização da distribuição espacial da variabilidade.
- **Análise em linha:** comparação dos valores estatísticos ao longo de uma linha horizontal da imagem.
- **Normalização:** comparação dos desvios padrão em relação a uma região de referência.
- **Filtro de valores extremos:** identificação e exclusão de regiões conforme critérios estatísticos configuráveis.
- **Exportação de resultados:** salvamento de mapas e gráficos em PNG e dos dados numéricos em CSV.
- **Interface gráfica:** configuração dos parâmetros e visualização dos resultados em um ambiente interativo.

## Metodologia

A análise é realizada por meio de uma janela móvel que percorre a imagem em posições sucessivas.

Para cada região de interesse, o programa calcula a média e o desvio padrão dos valores dos pixels. O desvio padrão é utilizado como métrica de variabilidade local, permitindo investigar como o ruído se distribui espacialmente.

Na análise bidimensional, os resultados são organizados em uma grade espacial e apresentados em um mapa de calor. Também é possível realizar uma análise ao longo de uma linha horizontal para observar a variação da métrica em função da posição.

No modo de mamografia, os resultados podem ser normalizados pelo desvio padrão de uma ROI de referência. Essa abordagem permite comparar a variabilidade relativa entre diferentes regiões da imagem.

O aplicativo também oferece um filtro baseado em percentis para identificar regiões com valores elevados de média ou desvio padrão. A influência desse filtro deve ser considerada na interpretação dos resultados.

## Tecnologias utilizadas

- **Python:** desenvolvimento da aplicação e implementação dos métodos de análise.
- **NumPy:** operações numéricas e cálculos estatísticos.
- **pydicom:** leitura de arquivos DICOM.
- **Matplotlib:** geração de gráficos e mapas de calor.
- **CustomTkinter:** construção da interface gráfica.
- **Pillow:** manipulação e apresentação de imagens.
- **CSV:** exportação de resultados para análise posterior.

## Como executar

### 1. Pré-requisitos

É necessário ter Python instalado e disponível no terminal.

### 2. Clonar o repositório

```bash
git clone https://github.com/lucasbarbosaas/mammostationarity.git
cd mammostationarity
```

### 3. Instalar as dependências

```bash
python -m pip install customtkinter pydicom numpy matplotlib pillow
```

### 4. Executar o aplicativo

```bash
python "Aplicativo estacionariedade final.py"
```

Após iniciar o programa, selecione uma imagem compatível, configure os parâmetros da análise e execute a varredura. Os resultados podem ser visualizados na interface e exportados pelos controles disponíveis.

## Resultados e interpretação

Os mapas de calor permitem observar diferenças espaciais nos valores de desvio padrão calculados nas ROIs. Os gráficos em linha complementam essa visualização, mostrando como a métrica varia ao longo de uma região específica da imagem.

A interpretação deve considerar as características da imagem, o tamanho das ROIs, o passo de varredura, a região de referência e os critérios de filtragem utilizados.

Uma variação espacial do desvio padrão, isoladamente, não determina a causa da variação nem comprova todas as propriedades estatísticas associadas à estacionariedade do ruído. Os resultados devem ser avaliados considerando o contexto físico e experimental.

## Objetivos do projeto

- Desenvolver uma ferramenta para análise quantitativa de imagens médicas.
- Automatizar cálculos estatísticos realizados em regiões de interesse.
- Facilitar a visualização e a comparação de padrões espaciais.
- Integrar programação, estatística e processamento de imagens em uma aplicação de pesquisa.
- Organizar os resultados em formatos que permitam análises complementares.

## Contexto acadêmico

O projeto está relacionado ao estudo da estacionariedade do ruído em sistemas digitais de imagens por raios X, no âmbito da Física Médica.

Seu desenvolvimento envolve a implementação de métodos numéricos, o processamento de imagens, a análise estatística e a visualização de dados para investigar características de sistemas de aquisição de imagens médicas.

## Limitações e uso responsável

O aplicativo tem finalidade de pesquisa e análise experimental. Seus resultados dependem dos dados de entrada, dos parâmetros selecionados e dos métodos empregados.

A ferramenta não substitui procedimentos de controle de qualidade validados nem deve ser utilizada isoladamente para diagnóstico ou decisões clínicas.

Ao trabalhar com imagens médicas, utilize dados autorizados para essa finalidade e remova informações pessoais identificáveis antes de compartilhar arquivos.

## Autor

**Lucas Barbosa Araújo da Silva**

Graduando em Física Médica na UFRJ, com interesse em programação, estatística, análise de dados e desenvolvimento de ferramentas computacionais.

[Repositório no GitHub](https://github.com/lucasbarbosaas/mammostationarity)
