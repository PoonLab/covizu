// utils.spec.js is created with Cypress, https://www.cypress.io/
// direct download of Cypress is required to run the tests

describe('isAccn Function', () => {
    it('Returns true if the string is an accession number', () => {
        cy.window().then((win) => {
            expect(win.isAccn("EPI_ISL_545540")).to.eq(true)
            expect(win.isAccn("EPI_ISL5455402")).to.eq(false)
        })
    })
})

describe('isLineage Function', () => {
    it('Retruns true if the string is a lineage ', () => {
        cy.window().then((win) => {
            expect(win.isLineage("B.101.511")).to.eq(false)
            expect(win.isLineage("B.1.511")).to.eq(true)
            expect(win.isLineage("B.1.585")).to.eq(true)
        })
    })
})

describe('formatDate Function', () => {
    it('Returns a string in an ISO8601 format', () => {
        cy.window().then((win) => {
            expect(win.formatDate(new Date('2021-07-05'))).to.eq('2021-07-05')
        })
    })
})

describe('isDate Function', () => {
    it('Returns true if the date is in the correct format (YYYY-MM-DD)', () => {
        cy.window().then((win) => {
            expect(win.isDate("Jul 05 2021")).to.eq(false)
            expect(win.isDate("2021-14-12")).to.eq(false)
            expect(win.isDate("2021-02-35")).to.eq(false)
            // expect(win.isDate("2021-02-31")).to.eq(false)
            expect(win.isDate("2021-07-05")).to.eq(true)
        })
    })
})