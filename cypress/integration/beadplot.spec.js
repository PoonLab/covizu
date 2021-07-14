// beadplot.spec.js is created with Cypress, https://www.cypress.io/
// direct download of Cypress is required to run the tests
import clusters from '../fixtures/test_clusters.json'

describe('Array-input functions', () => {
    let coldates = []
    it('Unique', () => {
        cy.window().then((win) => {
            Object.values(clusters[0]["nodes"]).forEach(e => e.forEach(e => coldates.push(e[0])))
            expect(win.unique(coldates).length).to.eq(12)
            expect(win.unique(coldates)).to.eql(['2020-05-02','2020-05-15', '2020-05-26', '2020-06-10', '2020-06-16', 
            '2020-06-05', '2020-07-05', '2020-07-07', '2020-06-03', '2020-06-24', '2020-07-01', '2020-07-08'])
        })
    })
    it('Mode', () => {
        cy.window().then((win) => {
            expect(win.mode(coldates)).to.eq("2020-05-15")
        })
    })
    it('Tabulate', () => {
        cy.window().then((win) => {
            expect(win.tabulate(coldates)).to.eql({'2020-05-02': 1,'2020-05-15': 2, '2020-05-26': 1, '2020-06-10': 1, '2020-06-16': 1, 
            '2020-06-05': 1, '2020-07-05': 1, '2020-07-07': 1, '2020-06-03': 1, '2020-06-24': 1, '2020-07-01': 1, '2020-07-08': 1})
        })
    })
})

describe('Parse Functions', () => {
    var result = [];
    var beaddata; 
    it('parse_variant', () => {
        cy.window().then((win) => {
            let variant = Object.values(clusters[0].nodes)[0];
            let accn = Object.keys(clusters[0].nodes)[0];
            result = win.parse_variant(variant, 0, 0, accn, null, null); 

            expect(result.variant).to.eql({
                'accession': "EPI_ISL_542678",
                'label': "USA/TX-HMH-1371",
                'x1': new Date('2020-05-02'),  
                'x2': new Date('2020-06-10'),  
                'y1': 0,
                'y2': 0,
                'count': 5,
                'country': {"USA": 5},
                'region': ['North America','North America','North America','North America','North America'],
                'numBeads': 4,
                'parent': null,
                'dist': 0,
                'unsampled': false
            })
            expect(result.points[1]).to.eql({
                'cidx': 0,
                'variant': "EPI_ISL_542678",
                'x': new Date('2020-05-15'),
                'y': 0,
                'count': 2,
                'accessions': ["EPI_ISL_545194", "EPI_ISL_545201"],
                'labels': ['USA/TX-HMH-2234', 'USA/TX-HMH-2246'],
                'region1': 'North America',
                'region': ['North America','North America'],
                'country': {"USA": 2},
                'parent': null,
                'dist': 0})
            })
    })
    it('parse_edgelist', () => {
        cy.window().then((win) => {
            // One variant 
            expect(win.parse_edgelist(clusters[0], [result.variant], [result.points])).to.eql([])

            // More than one variant
            let variants = [result.variant]
            variants.push({
                'accession': "EPI_ISL_543301",
                'label': "TX-HMH-4019",
                'x1': new Date('2020-06-16'),  
                'x2': new Date('2020-06-16'),  
                'y1': 0,
                'y2': 0,
                'count': 1,
                'country': {"USA": 1},
                'region': ['North America'],
                'numBeads': 1,
                'parent': null,
                'dist': 0,
                'unsampled': false
            })

            let points = [result.points]
            points = points.concat({
                'cidx': 1,
                'variant': "EPI_ISL_543301",
                'x': new Date("2020-06-16"),
                'y': 0,
                'count': 1,
                'accessions': ["EPI_ISL_545450"],
                'labels': ['USA/TX-HMH-2649'],
                'region1': 'North America',
                'region': ['North America'],
                'country': {"USA": 1},
                'parent': null,
                'dist': 0
            })

            expect(win.parse_edgelist(clusters[0], variants, points)).to.eql([
                {
                    'y1': 0,
                    'y2': 0,
                    'x1': new Date("2020-06-16"),  
                    'x2': new Date("2020-06-16"),
                    'parent': "USA/TX-HMH-1371",
                    'child': "TX-HMH-4019",
                    'dist': 3.12,
                    'support': undefined
                }
            ])
        })
    })
    it('parse_clusters', () => {
        // Passing of the test is dependent on parse_edgelist and parse_variant functions
        cy.window().then((win) => {
            let edgelist = [{'y1': 1, 'y2': 2, 'x1': new Date("2020-06-16"), 'x2': new Date("2020-06-16"), 'parent': "USA/TX-HMH-1371", 'child': "USA/TX-HMH-4019", 'dist': 3.12, 'support': undefined},
            {'y1': 1, 'y2': 3, 'x1': new Date("2020-05-02"), 'x2': new Date("2020-05-02"), 'parent': "USA/TX-HMH-1371", 'child': "unsampled0", 'dist': 1.61, 'support': 0.62},
            {'y1': 1, 'y2': 4, 'x1': new Date("2020-06-05"), 'x2': new Date("2020-06-05"), 'parent': "USA/TX-HMH-1371", 'child': "USA/TX-HMH-2649", 'dist': 2.67, 'support': undefined},
            {'y1': 1, 'y2': 5, 'x1': new Date("2020-07-05"), 'x2': new Date("2020-07-05"), 'parent': "USA/TX-HMH-1371", 'child': "USA/TX-HMH-MCoV-10847", 'dist': 3.36, 'support': 0.94},
            {'y1': 3, 'y2': 6, 'x1': new Date("2020-06-03"), 'x2': new Date("2020-06-03"), 'parent': "unsampled0", 'child': "USA/TX-HMH-2767", 'dist': 1.24, 'support': undefined},
            {'y1': 3, 'y2': 7, 'x1': new Date("2020-06-24"), 'x2': new Date("2020-06-24"), 'parent': "unsampled0", 'child': "USA/TX-DSHS-1864", 'dist': 2.72, 'support': undefined},
            {'y1': 5, 'y2': 8, 'x1': new Date("2020-07-01"), 'x2': new Date("2020-07-01"), 'parent': "USA/TX-HMH-MCoV-10847", 'child': "USA/TX-HMH-MCoV-8639", 'dist': 0.87, 'support': undefined},
            {'y1': 5, 'y2': 9, 'x1': new Date("2020-07-08"), 'x2': new Date("2020-07-08"), 'parent': "USA/TX-HMH-MCoV-10847", 'child': "USA/TX-HMH-MCoV-10410", 'dist': 1.12, 'support': undefined}]

            let variants = [{"accession": "EPI_ISL_542678","label": "USA/TX-HMH-1371","x1": new Date("2020-05-02"),"x2": new Date("2020-07-05"),"y1": 1,"y2": 1,"count": 5,"country": {"USA": 5},
                    "region": [ "North America", "North America", "North America", "North America", "North America"],"numBeads": 4,"parent": null,"dist": 0,"unsampled": false },
                {"accession": "EPI_ISL_543301","label": "USA/TX-HMH-4019","x1": new Date("2020-06-16"),"x2": new Date("2020-06-16"),"y1": 2,"y2": 2,"count": 1,"country": {"USA": 1},
                    "region": [ "North America"],"numBeads": 1,"parent": "USA/TX-HMH-1371","dist": 3.12,"unsampled": false},
                {"accession": "unsampled0","label": "unsampled0","x1": new Date("2020-05-02"),"x2": new Date("2020-07-08"),"y1": 3,"y2": 3,"count": 0,"country": null,
                    "region": null,"numBeads": 0,"parent": "USA/TX-HMH-1371","dist": 1.61,"unsampled": true},
                {"accession": "EPI_ISL_545450","label": "USA/TX-HMH-2649","x1": new Date("2020-06-05"),"x2": new Date("2020-06-05"),"y1": 4,"y2": 4,"count": 1,"country": {"USA": 1},
                    "region": ["North America"],"numBeads": 1,"parent": "USA/TX-HMH-1371","dist": 2.67,"unsampled": false },
                {"accession": "EPI_ISL_789657","label": "USA/TX-HMH-MCoV-10847","x1": new Date("2020-07-01"),"x2": new Date("2020-07-08"),"y1": 5,"y2": 5,"count": 2,"country": {"USA": 2},
                    "region": ["North America","North America"],"numBeads": 2,"parent": "USA/TX-HMH-1371","dist": 3.36,"unsampled": false},
                {"accession": "EPI_ISL_545540","label": "USA/TX-HMH-2767","x1": new Date("2020-06-03"),"x2": new Date("2020-06-03"),"y1": 6,"y2": 6,"count": 1,"country": {"USA": 1},
                    "region": ["North America"],"numBeads": 1,"parent": "unsampled0","dist": 1.24,"unsampled": false},
                {"accession": "EPI_ISL_745299","label": "USA/TX-DSHS-1864","x1": new Date("2020-06-24"),"x2": new Date("2020-06-24"),"y1": 7,"y2": 7,"count": 1,"country": {    "USA": 1},
                    "region": ["North America"],"numBeads": 1,"parent": "unsampled0","dist": 2.72,"unsampled": false},
                {"accession": "EPI_ISL_780668","label": "USA/TX-HMH-MCoV-8639","x1": new Date("2020-07-01"),"x2": new Date("2020-07-01"),"y1": 8,"y2": 8,"count": 1,"country": {    "USA": 1},
                    "region": ["North America"],"numBeads": 1,"parent": "USA/TX-HMH-MCoV-10847","dist": 0.87,"unsampled": false},
                {"accession": "EPI_ISL_789458","label": "USA/TX-HMH-MCoV-10410","x1": new Date("2020-07-08"),"x2": new Date("2020-07-08"),"y1": 9,"y2": 9,"count": 1,"country": {"USA": 1},
                    "region": ["North America"],"numBeads": 1,"parent": "USA/TX-HMH-MCoV-10847","dist": 1.12,"unsampled": false}
            ]

            let points = [{"cidx": "0","variant": "EPI_ISL_542678","x": new Date("2020-05-02"),"y": 1,"count": 1,"accessions": [
                        "EPI_ISL_542678"],"labels": ["USA/TX-HMH-1371"],"region1": "North America","region": ["North America"],"country": {"USA": 1},
                        "parent": null,"dist": 0},
                        {"cidx": "0","variant": "EPI_ISL_542678","x": new Date("2020-05-15"),"y": 1,"count": 2,"accessions": 
                        ["EPI_ISL_545194","EPI_ISL_545201"],"labels": ["USA/TX-HMH-2234","USA/TX-HMH-2246"],"region1": "North America","region": ["North America","North America"],
                        "country": {"USA": 2},"parent": null,"dist": 0},
                        {"cidx": "0","variant": "EPI_ISL_542678","x": new Date("2020-05-26"),"y": 1,"count": 1,"accessions": 
                        ["EPI_ISL_545290"],"labels": ["USA/TX-HMH-2403"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": null,"dist": 0},
                        {"cidx": "0","variant": "EPI_ISL_542678","x": new Date("2020-06-10"),"y": 1,"count": 1,"accessions": 
                        ["EPI_ISL_543906"],"labels": ["USA/TX-HMH-2976"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1}, "parent": null,"dist": 0},
                        {"cidx": "0","variant": "EPI_ISL_543301","x": new Date("2020-06-16"),"y": 2,"count": 1,"accessions": 
                        ["EPI_ISL_543301"],"labels": ["USA/TX-HMH-4019"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": "USA/TX-HMH-1371","dist": 3.12},
                        {"cidx": "0","variant": "EPI_ISL_545450","x": new Date("2020-06-05"),"y": 4,"count": 1,"accessions": 
                        ["EPI_ISL_545450"],"labels": ["USA/TX-HMH-2649"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": "USA/TX-HMH-1371","dist": 2.67},
                        {"cidx": "0","variant": "EPI_ISL_789657","x": new Date("2020-07-05"),"y": 5,"count": 1,"accessions": 
                        ["EPI_ISL_789657"],"labels": ["USA/TX-HMH-MCoV-10847"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": "USA/TX-HMH-1371","dist": 3.36},
                        {"cidx": "0","variant": "EPI_ISL_789657","x": new Date("2020-07-07"),"y": 5,"count": 1,"accessions": 
                        ["EPI_ISL_784276"],"labels": ["USA/TX-HMH-MCoV-8225"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": "USA/TX-HMH-1371","dist": 3.36},
                        {"cidx": "0","variant": "EPI_ISL_545540","x": new Date("2020-06-03"),"y": 6,"count": 1,"accessions": 
                        ["EPI_ISL_545540"],"labels": ["USA/TX-HMH-2767"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": "unsampled0","dist": 1.24},
                        {"cidx": "0","variant": "EPI_ISL_745299","x": new Date("2020-06-24"),"y": 7,"count": 1,"accessions": 
                        ["EPI_ISL_745299"],"labels": ["USA/TX-DSHS-1864"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": "unsampled0","dist": 2.72},
                        {"cidx": "0","variant": "EPI_ISL_780668","x": new Date("2020-07-01"),"y": 8,"count": 1,"accessions": 
                        ["EPI_ISL_780668"],"labels": ["USA/TX-HMH-MCoV-8639"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": "USA/TX-HMH-MCoV-10847","dist": 0.87},
                        {"cidx": "0","variant": "EPI_ISL_789458","x": new Date("2020-07-08"),"y": 9,"count": 1,"accessions": 
                        ["EPI_ISL_789458"],"labels": ["USA/TX-HMH-MCoV-10410"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": "USA/TX-HMH-MCoV-10847","dist": 1.12}]

            beaddata = win.parse_clusters(clusters)
            expect(beaddata[0].edgelist).to.eql(edgelist)
            expect(beaddata[0].points).to.eql(points)
            expect(beaddata[0].variants).to.eql(variants)
                
        })
    })
    it('Merge tables ', () => {
        cy.window().then((win) => {
            expect(win.merge_tables(beaddata[0].variants.map(x => x.country))).to.eql({'USA': 13})
        })
    })
    it('Beadplot function', () => {
        cy.window().then((win) => {
            // Checks that beadplot corresponds to selected lineage
            win.beadplot(786)
            cy.get('#svg-clusteraxis>svg>g').then(($axis) => {
                cy.get('.beadplot-content>svg>g').then(($beadplot) => {
                    cy.get('[cidx="cidx-786"]').click()
                    cy.get($axis).children().each(($el) => {
                        cy.get('#svg-clusteraxis>svg>g').find($el)
                    })
                    cy.get($beadplot).children().each(($el) => {
                        cy.get('.beadplot-content>svg>g').find($el)
                    })
                })
            })
        })
    })
})

describe('Tables', () => {
    it('region_to_string', () => {
        cy.window().then((win) => {
            expect(win.region_to_string({'North America':5,'Asia':10})).to.eq("<b>Number of cases:</b><br>&nbsp;&nbsp;North America: 5<br>&nbsp;&nbsp;Asia: 10<br>Total: 15<br>")
        })
    })
    it('Sequence details: gen_details_table', () => {
        cy.contains('B.1.511').click()
        cy.get('circle:visible').first().click().invoke('attr','id').as('lineage_id')
        cy.get('@lineage_id').then(id => {
            cy.get('.selectionH').should('have.attr', 'bead', `${id}`)
        })

        cy.get('#tabs').contains('Samples').click()
        cy.get('#tabs-2').contains('2020-05-01')
        cy.get('#tabs-2').contains('USA/TX-HMH-1371')
        cy.get('#tabs-2').contains('EPI_ISL_542678')
    }) 
})